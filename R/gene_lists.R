library(data.table)
library(randomForest)
library(MiRaGE)
library(Biobase)
library(org.Hs.eg.db)
library(SummarizedExperiment)
library(DT)

setwd("C:/Users/Dovydas/Desktop/gbm/R/gene_analysis/Gene Analysis/")
source("./functions.R")
hg_file <- "./data/HG-U133A_2018-02-09_TCGA_GBM_expression.txt"
agilent_file <- "./data/Agilent-4502A_2018-02-09_TCGA_GBM_expression.txt"
pheno_file <- "./data/2018-02-09_TCGA_GBM_pheno.txt"  

alias_to_eg <- as.list(org.Hs.egALIAS2EG)
refseq <- as.list(org.Hs.egREFSEQ)

#hg <- readRDS(file="./HG-U133A/HG-U133A_data_complete.rds")
#agilent <- readRDS(file="./Agilent-4502A/Agilent-4502A_data_complete.rds")
glist <- scan(file="./gene_overlap_list.txt", what=character())

hg <- read_exp_data(pheno_file, hg_file, "Subtype_Verhaak_2010")
#agilent <- read_exp_data(pheno_file, agilent_file, "Subtype_Verhaak_2010")

hg_red <- hg[, c(Reduce(intersect, list(glist, colnames(hg))), "Subtype")]
#agilent_red <- agilent[, c(Reduce(intersect, list(glist, colnames(agilent))), "Subtype")]

clinical_data <- read.csv(pheno_file, sep="\t")

hg_control <- hg_red[is.na(hg_red$Subtype), ]
hg_control$Subtype <- NULL

#hg_red <- hg_red[rownames(hg_red) %in% treatment_types$Sample, ]

#agilent_control <- agilent_red[is.na(agilent_red$Subtype), ]
#agilent_control$Subtype <- NULL

#Reduce(intersect, list(glist, colnames(agilent)))
#write.csv(hg_red, file="./hg.csv")
#write.csv(agilent_red, file="./agilent.csv")

# Random Forest -----------------------------------------------------------

train.prop <- 0.9
data.complete <- hg_red
data.classical <- subset(data.complete, Subtype=="Classical")
data.classical <- data.classical[sample(nrow(data.classical)), ]
data.mesenchymal <- subset(data.complete, Subtype=="Mesenchymal")
data.mesenchymal <- data.mesenchymal[sample(nrow(data.mesenchymal)), ]
data.neural <- subset(data.complete, Subtype=="Neural")
data.neural <- data.neural[sample(nrow(data.neural)), ]
data.proneural <- subset(data.complete, Subtype=="Proneural")
data.proneural <- data.proneural[sample(nrow(data.proneural)), ]

train <- data.classical[1:(train.prop * nrow(data.classical)),]
train <- rbind(train, data.mesenchymal[1:(train.prop * nrow(data.mesenchymal)),])
train <- rbind(train, data.neural[1:(train.prop * nrow(data.neural)),])
train <- rbind(train, data.proneural[1:(train.prop * nrow(data.proneural)),])
train <- train[sample(nrow(train)), ]

test <- data.classical[(train.prop * nrow(data.classical)):nrow(data.classical),]
test <- rbind(test, data.mesenchymal[(train.prop * nrow(data.mesenchymal)):nrow(data.mesenchymal),])
test <- rbind(test, data.neural[(train.prop * nrow(data.neural)):nrow(data.neural),])
test <- rbind(test, data.proneural[(train.prop * nrow(data.proneural)):nrow(data.proneural),])
test <- test[sample(nrow(test)), ]

model_rf <- randomForest(Subtype ~ ., data = train)
preds <- predict(model_rf, test)
table(preds)

sum((unname(preds) == test$Subtype) == TRUE)/length(test$Subtype)

# MiRaGE ------------------------------------------------------------------

treatment_types <- clinical_data[grepl("Alkylating", clinical_data$Therapy_Class, fixed = TRUE), c("Sample", "Therapy_Class")]
treatment_types <- clinical_data[clinical_data$Therapy %in% c("TMZ Chemoradiation, TMZ Chemo"), c("Sample", "Therapy_Class")]
exp_data <- data.classical
exp_data <- hg_red[rownames(hg_red) %in% treatment_types$Sample, ]
control_data <- data.mesenchymal[, 1:(ncol(data.mesenchymal)-1)]
#subtype_data <- data.mesenchymal

genes <- colnames(exp_data)[1:(ncol(exp_data)-1)]

eg <- alias_to_eg[genes]

genexp <- exp_data[, 1:(ncol(exp_data)-1), drop=F]
genexp <- rbind(control_data, genexp)
genexp <- as.data.frame(t(genexp))
rownames(genexp) <- c()
genexp <- cbind(genes, genexp)

control_labels <- c()
for(i in c(1:nrow(control_data))){
  control_labels <- c(control_labels, paste("con.", i, sep=""))
}

labels <- c()
for(i in c(1:nrow(exp_data))){
  labels <- c(labels, paste("gbm.", i, sep=""))
}

names(genexp) <- c("gene", control_labels, labels)

new_df <- data.frame(Gene=character(), EG=character())
pos <- 1
for(i in unname(eg)){
  gen <- names(eg)[pos]
  for(j in i){
    new_df <- rbind(new_df, data.frame(Gene=gen, 
                                       EG=j, 
                                       genexp[genexp$gene==gen, 2:ncol(genexp)]))
  }
  pos <- pos + 1
}

eg_list <- new_df$EG
ref <- refseq[eg_list]

df <- data.frame(RefSeq = character())
pos <- 1
for(i in unname(ref)){
  for(j in i){
    if(j %like% "NM"){
      gen <- as.character(new_df$Gene[pos])
      df <- rbind(df, data.frame(RefSeq=j,
                                 genexp[genexp$gene == gen, 2:ncol(genexp)]))
    }
  }
  pos <- pos + 1
}

rownames(df) <- c(1:nrow(df))

expset <- new("ExpressionSet", expr=as.matrix(df[,-1]))
fData(expset)[["gene_id"]] <- df[,1]
pData(expset)[["sample_name"]] <- colnames(df)[-1]

results <- MiRaGE(expset, species="HS")

results$P0[order(results$P0[,2])[1:20], ]
#results$P1[order(results$P1[,2])[1:20], ]

# Result Writing ----------------------------------------------------------

write.csv(results$P0[order(results$P0[,2])[1:20], ], "./MiRaGE Results/tmz_chemoradiation_tmz_chemo_up.csv", row.names=F)
write.csv(results$P1[order(results$P1[,2])[1:20], ], "./MiRaGE Results/tmz_chemoradiation_tmz_chemo_down.csv", row.names=F)

gbm_up <- read.csv("./MiRaGE Results/HG-U133A/gbm_up.csv")
classical_up <- read.csv("./MiRaGE Results/HG-U133A/classical_up.csv")
mesenchymal_up <- read.csv("./MiRaGE Results/HG-U133A/mesenchymal_up.csv")
neural_up <- read.csv("./MiRaGE Results/HG-U133A/neural_up.csv")
proneural_up <- read.csv("./MiRaGE Results/HG-U133A/proneural_up.csv")

gbm_vs_classical <- Reduce(setdiff, list(gbm_up$Refseq, classical_up$Refseq))
gbm_vs_mesenchymal <- Reduce(setdiff, list(gbm_up$Refseq, mesenchymal_up$Refseq))

Reduce(setdiff, list(gbm_up$Refseq, neural_up$Refseq))
Reduce(setdiff, list(gbm_up$Refseq, proneural_up$Refseq))

Reduce(setdiff, list(gbm_vs_classical, gbm_vs_mesenchymal))
