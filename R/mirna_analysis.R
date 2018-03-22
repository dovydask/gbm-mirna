library(data.table)
library(combinat)
library(ggbiplot)
library(rpart)
library(pracma)
library(factoextra)
library(FactoMineR)
library(party)
library(relaimpo)
library(randomForest)
library(Metrics)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(reshape)
library(scales)

setwd("C:/Users/Dovydas/Desktop/gbm/R/mirna_analysis")
pheno.file <- "./data/2018-02-09_TCGA_GBM_pheno.txt"
mirna.files <- Sys.glob("./data/GBM - miRNA Expression/*.txt")
control.files <- Sys.glob("./data/Control-patients/Expression-miRNA/UNC__H-miRNA_8x15K/Level_3/*.txt")

# Data Frame Generation ---------------------------------------------------

first_file <- read.delim(mirna.files[1])
mirna.list <- first_file$miRNA.id
mirna.data <- data.frame(matrix(ncol=(length(mirna.list)+2), nrow=0))
colnames(mirna.data) <- c("Case", "Sample", as.character(mirna.list))
for (file in mirna.files){
  print(file)
  case.data <- read.delim(file)
  case <- case.data$barcode[1]
  sample <- gsub("-", ".", substr(case.data$barcode[1], 0, 12))
  stopifnot(case.data$miRNA.id == mirna.list)
  mirna.data[nrow(mirna.data) + 1, ] = c(as.character(case), sample, as.numeric(case.data$value))
}

pheno.data <- read.delim(pheno.file)
mirna.data <- mirna.data[mirna.data$Sample %in% pheno.data$Sample, ]
mirna.data["Subtype"] <- pheno.data$Subtype_Verhaak_2010[match(mirna.data$Sample, pheno.data$Sample)]
mirna.data <- mirna.data[complete.cases(mirna.data),]
rownames(mirna.data) <- mirna.data$Case
mirna.data$Case <- NULL
write.csv(mirna.data, file="./mirna_data.csv")
mirna.data.original <- mirna.data

# Control Patients --------------------------------------------------------

control.first_file <- read.delim(control.files[1])
control.mirna.list <- control.first_file$miRNA.id
control.mirna.data <- data.frame(matrix(ncol=(length(control.mirna.list)+2), nrow=0))
colnames(control.mirna.data) <- c("Case", "Sample", as.character(control.mirna.list))
for (file in control.files){
  print(file)
  case.data <- read.delim(file)
  case <- case.data$barcode[1]
  sample <- gsub("-", ".", substr(case.data$barcode[1], 0, 12))
  stopifnot(case.data$miRNA.id == control.mirna.list)
  control.mirna.data[nrow(control.mirna.data) + 1, ] = c(as.character(case), sample, as.numeric(case.data$value))
}

# Treatment ---------------------------------------------------------------

rad_treatment <- pheno.data[pheno.data$Therapy_Class %in% c("Standard Radiation",
                                                          "Nonstandard Radiation",
                                                          "Unspecified Radiation"),
                            c("Sample", "CIMP_status", "IDH1_status", "MGMT_status", "Therapy_Class")]


rad_gbm <- mirna.data.original[mirna.data.original$Sample %in% rad_treatment$Sample, ]
rad_gcimp <- rad_treatment[rad_treatment$CIMP_status == "G-CIMP", ]
rad_gcimp <- mirna.data.original[mirna.data.original$Sample %in% rad_gcimp$Sample, ]


# Linear Regression -------------------------------------------------------

dict <- vector(mode="list", length=4)
names(dict) <- c("Classical", "Mesenchymal", "Neural", "Proneural")
dict[[1]] <- 1
dict[[2]] <- 2
dict[[3]] <- 3
dict[[4]] <- 4

mirna.data <- rad_gcimp
max_pvalue <- 0.025
mirna.lr <- data.table(miRNA=character(), pValue=numeric())
i <- 1
eh <- lapply(levels(factor(mirna.data$Subtype)), function(x) dict[x])
class_levels <- unname(unlist(eh))
classes <- as.numeric(mirna.data$Subtype)
for (miRNA in names(mirna.data)[2:(length(names(mirna.data))-1)]) {
  print(i)
  i <- i + 1
  if (miRNA != "Subtype") {
    for (permutation in permn(class_levels)) {
      mirna.expressions <- as.numeric(mirna.data[[miRNA]])
      df <- data.frame(classes, mirna.expressions)
      df$classes <- as.numeric(factor(df$classes, levels=permutation))
      df <- df[order(df$classes), ]
      model <- lm(df$mirna.expressions ~ df$classes)
      pvalue <- as.numeric(anova(model)$'Pr(>F)'[1])
      if (pvalue < max_pvalue) {
        mirna.lr <- rbindlist(list(mirna.lr, list(miRNA, pvalue)))
      }
    }
  }
}

list_length <- 150
pvalue_cutoff <- 0.025

colnames(mirna.lr) <- c("miRNA", "pValue")
mirna.ordered <- mirna.lr[order(mirna.lr$miRNA, mirna.lr$pValue), ]
mirna.unique <- subset(mirna.ordered, !duplicated(miRNA))
mirna.filtered <- mirna.unique[mirna.unique$pValue < pvalue_cutoff, ]
mirna.filtered <- mirna.filtered[order(mirna.filtered$pValue), ]
mirna.subtyped <- subset(mirna.data, select=c(as.vector(mirna.filtered$miRNA), "Subtype"))
mirna.complete <- mirna.subtyped[complete.cases(mirna.subtyped), ]
#write(mirna.filtered$miRNA[1:list_length], file="./mirna_lm_list.txt")

# Linear Discriminant Analysis --------------------------------------------

interest_list <- scan("./mirna_list.txt", what=character())
interest_data <- mirna.data[, c(interest_list, "Subtype")]

#lab_list <- c("hsa-miR-181d", "hsa-miR-21", "hsa-miR-30b", "hsa-miR-30c", "hsa-miR-425-5p", "hsa-miR-93",
#              "hsa-miR-10b", "hsa-miR-124a", "hsa-miR-129", "hsa-miR-139", "hsa-miR-15b", "hsa-miR-218", 
#              "hsa-miR-7", "hsa-miR-34a", "hsa-miR-155")
#interest_data <- mirna.data.original[, c(lab_list, "Subtype")]

colnames(interest_data) <- make.names(colnames(interest_data), unique = T)
cols <- c(1:(ncol(interest_data)-1))
interest_data[,cols] <- apply(interest_data[,cols], 2, function(x) as.numeric(as.character(x)));

interest_data_scaled <- as.data.frame(scale(interest_data[, 1:(ncol(interest_data)-1)]))
interest_data_scaled["Subtype"] <- interest_data["Subtype"]

train.prop <- 0.9
data.complete <- interest_data_scaled
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

#train <- data.neural[1:(train.prop * nrow(data.neural)), ]
#train <- rbind(train, data.proneural[1:(train.prop * nrow(data.proneural)),])
#train$Subtype <- factor(train$Subtype)

#test <- data.neural[(train.prop * nrow(data.neural)):nrow(data.neural), ]
#test <- rbind(train, data.proneural[(train.prop * nrow(data.proneural)):nrow(data.proneural),])
#test$Subtype <- factor(test$Subtype)

model_rf <- randomForest(Subtype ~ ., data = train)
preds <- predict(model_rf, test)
table(preds)

sum((unname(preds) == test$Subtype) == TRUE)/length(test$Subtype)

importance <- as.data.frame(importance(model_rf))
importance <- importance[order(-importance$MeanDecreaseGini), , drop=FALSE]

top_list_count <- 50
top <- paste(rownames(importance)[1:top_list_count], collapse = '+')
formula <- paste("Subtype ~ ", top)
model_rf_imp <- randomForest(as.formula(formula), data=train)
preds_imp <- predict(model_rf_imp, test)
sum((unname(preds_imp) == test$Subtype) == TRUE)/length(test$Subtype)

#write(rownames(importance)[1:top_list_count], file="./best_predictors.txt")

top_data <- interest_data[, c(rownames(importance)[1:top_list_count], "Subtype")]
top_data <- top_data[order(top_data$Subtype), ]
plot_data <- scale(as.matrix(top_data[,1:(ncol(top_data)-1)]))

# Heatmap -----------------------------------------------------------------

library(ComplexHeatmap)
library(circlize)
plot_data_t <- t(plot_data)
ha <- HeatmapAnnotation(df = data.frame(type = top_data$Subtype),
                        col = list(type = c("Classical" = "red", "Mesenchymal" = "forestgreen",
                                            "Neural" = "blue", "Proneural" = "yellow"))
                        )

Heatmap(plot_data_t, cluster_columns = F, 
        show_column_dend = FALSE, show_column_names = FALSE,
        km = 4, top_annotation = ha
        )
# Plots -------------------------------------------------------------------

file.remove(file.path("./plots/all/", list.files("./plots/all/")))
file.remove(file.path("./plots/rejected/", list.files("./plots/rejected/")))
file.remove(file.path("./plots/confirmed/", list.files("./plots/confirmed/")))
file.remove(file.path("./plots/control/", list.files("./plots/control")))

classes <- mirna.data$Subtype
for (mirna in colnames(mirna.data)[3:(ncol(mirna.data)-1)]){
  mirna.expressions <- as.numeric(mirna.data[[mirna]])
  png_filename <- paste("./plots/all/", mirna, ".png", sep="")
  png_filename <- gsub("\\*", "[a]", png_filename)
  png(png_filename)
  plot(classes, mirna.expressions, ylim = c(5, 15))
  dev.off()
}

classes <- mirna.rejected$Subtype
mirna.rejected <- mirna.data[, !(names(mirna.data) %in% colnames(mirna.complete)[1:(length(colnames(mirna.complete))-1)])]
for (mirna in colnames(mirna.rejected)[3:(ncol(mirna.rejected)-1)]){
  mirna.expressions <- as.numeric(mirna.rejected[[mirna]])
  png_filename <- paste("./plots/rejected/", mirna, ".png", sep="")
  png_filename <- gsub("\\*", "[a]", png_filename)
  png(png_filename)
  plot(classes, mirna.expressions, ylim = c(5, 15))
  dev.off()
}

classes <- mirna.complete$Subtype
for (mirna in colnames(mirna.complete)[3:(ncol(mirna.complete)-1)]){
  mirna.expressions <- as.numeric(mirna.complete[[mirna]])
  png_filename <- paste("./plots/confirmed/", mirna, ".png", sep="")
  png_filename <- gsub("\\*", "[a]", png_filename)
  png(png_filename)
  plot(classes, mirna.expressions, ylim = c(5, 15))
  dev.off()
}

for (mirna in colnames(control.mirna.data)[3:(ncol(control.mirna.data)-1)]){
  mirna.expressions <- as.numeric(control.mirna.data[[mirna]])
  png_filename <- paste("./plots/control/", mirna, ".png", sep="")
  png_filename <- gsub("\\*", "[a]", png_filename)
  png(png_filename)
  boxplot(mirna.expressions, ylim = c(5, 15))
  dev.off()
}



