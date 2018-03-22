#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")
library("TCGAWorkflow")
library("TCGAbiolinks")
library("TCGAWorkflowData")
library("SummarizedExperiment")
library("DT")
library("MiRaGE")
library("org.Hs.eg.db")
library("data.table")
library("plyr")

setwd("C:/Users/Dovydas/Desktop/gbm/R/mirna_analysis/TCGAWorkflow")

query.exp.gbm <- GDCquery(project = "TCGA-GBM",
                          legacy = TRUE,
                          data.category = "Gene expression",
                          data.type = "Gene expression quantification",
                          platform = "Illumina HiSeq",
                          file.type = "results",
                          sample.type = "Primary solid Tumor"
                          )
GDCdownload(query.exp.gbm, method="client")
exp.gbm <- GDCprepare(query = query.exp.gbm, 
                      save = TRUE, 
                      save.filename = "gbmExp.rda"
)

query.exp.gbm.control <- GDCquery(project = "TCGA-GBM",
                          legacy = TRUE,
                          data.category = "Gene expression",
                          data.type = "Gene expression quantification",
                          platform = "Illumina HiSeq",
                          file.type = "results",
                          sample.type = "Solid Tissue Normal"
                          )
GDCdownload(query.exp.gbm.control, method="client")
exp.gbm.control <- GDCprepare(query = query.exp.gbm.control, 
                              save = TRUE, 
                              save.filename = "gbmExpControl.rda",
                              summarizedExperiment = FALSE
                              )

exp.gbm.control <- exp.gbm.control[, c(1, 4, 7, 10, 13)]

data <- assay(exp.gbm)
genes.info <- rowRanges(exp.gbm)
sample.info <- colData(exp.gbm)

# datatable(as.data.frame(sample.info),
#           options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
#           rownames = FALSE
#           )

gbm.subtypes <- TCGAquery_subtype(tumor = "gbm")
gbm_clin <- GDCquery_clinic(project = "TCGA-GBM", type = "Clinical")

query.clinical.gbm <- GDCquery(project = "TCGA-GBM",
                               data.category = "Clinical")

GDCdownload(query.clinical.gbm, method = "client")
gbm_treatment_radiation <- GDCprepare_clinic(query.clinical.gbm, "radiation")
gbm_treatment_drug <- GDCprepare_clinic(query.clinical.gbm, "drug")
write.table(gbm_treatment_radiation, "./gbm_clin_radiation.txt", sep="\t", row.names=F)
write.table(gbm_treatment_drug, "./gbm_clin_drug.txt", sep="\t", row.names=F)

# datatable(gbm_clin[1:10,], 
#           options = list(scrollX = TRUE, keys = TRUE), 
#           rownames = FALSE
#           )
# datatable(gbm.subtypes[1:10,], 
#           options = list(scrollX = TRUE, keys = TRUE), 
#           rownames = FALSE
#           )

expression_data <- readRDS("./expression_data.rds")
expression_data <- assays(exp.gbm)$raw_count
expression_data <- unique(expression_data)

gbm.subtypes_df <- data.frame(patient = gbm.subtypes$patient, Original.Subtype = gbm.subtypes$Original.Subtype)
gbm.subtypes_df <- gbm.subtypes_df[gbm.subtypes_df$Original.Subtype %in% c("Classical", "Mesenchymal", "Neural", "Proneural"), ]
gbm.subtypes_df$Original.Subtype <- factor(gbm.subtypes_df$Original.Subtype, levels=c("Classical", "Mesenchymal", "Neural", "Proneural"))

new_names <- c()
for(i in rownames(expression_data)){
  new_names <- c(new_names, strsplit(i, split='|', fixed=TRUE)[[1]][1])
}

colnames(expression_data) <- sapply(colnames(expression_data), substring, 0, 12, USE.NAMES=F)

expression_data <- cbind(new_names, expression_data)
rownames(expression_data) <- c(1:nrow(expression_data))
expression_data <- as.data.frame(expression_data)

gbm.subtypes_df <- gbm.subtypes_df[gbm.subtypes_df$patient %in% colnames(expression_data), ]

saveRDS(gbm.subtypes_df, "./gbm_subtypes_df.rds")
gene_list <- scan("../../gene_analysis/Gene Analysis/TCGA/final_list.txt", what=character())
gene_list <- gene_list[2:length(gene_list)]

expression_data <- expression_data[, c("new_names", as.character(gbm.subtypes_df$patient))]
expression_data <- expression_data[order(expression_data$new_names), ]
expression_data <- expression_data[expression_data$new_names != "?", ]
expression_data <- expression_data[!duplicated(expression_data$new_names), ]
expression_data <- expression_data[expression_data$new_names %in% gene_list, ]
rownames(expression_data) <- c(1:nrow(expression_data))

control_data <- readRDS("./control_data.rds")
control_data <- exp.gbm.control
control_data <- unique(control_data)

new_control_names <- c()
for(i in rownames(control_data)){
  new_control_names <- c(new_control_names, strsplit(i, split='|', fixed=TRUE)[[1]][1])
}

control_data <- cbind(new_control_names, control_data)
rownames(control_data) <- c(1:nrow(control_data))
control_data <- control_data[order(control_data$new_control_names), ]
control_data <- control_data[control_data$new_control_names != "?", ]
control_data <- control_data[!duplicated(control_data$new_control_names), ]
control_data <- control_data[control_data$new_control_names %in% gene_list, ]
control_data[control_data$new_control_names %in% gene_list, ]

intersection_list <- Reduce(intersect, list(control_data$new_control_names, expression_data$new_names))
rownames(expression_data) <- expression_data$new_names
expression_data <- expression_data[c(intersection_list), ]
rownames(expression_data) <- c(1:nrow(expression_data))
rownames(control_data) <- control_data$new_control_names
control_data <- control_data[c(intersection_list), ]
rownames(control_data) <- c(1:nrow(control_data))

saveRDS(expression_data, file="./expression_data_reduced.rds")
saveRDS(control_data, file="./control_data_reduced.rds")

#expression_data <- readRDS("./expression_data.rds")
#control_data <- readRDS("./control_data.rds")

gbm_dataset <- cbind(control_data, expression_data[, 2:ncol(expression_data)])
colnames(gbm_dataset)[1:6] <- c("Gene", sapply(colnames(gbm_dataset)[2:6], substr, 11, 22, USE.NAMES = F))

alias_to_eg <- as.list(org.Hs.egALIAS2EG)
refseq <- as.list(org.Hs.egREFSEQ)
eg <- alias_to_eg[as.character(gbm_dataset$Gene)]

time <- Sys.time()
new_df <- data.frame(Gene=character(), EG=character())
pos <- 1
for(i in unname(eg)){
  gen <- names(eg)[pos]
  for(j in i){
    new_df <- rbind.fill(new_df, data.frame(Gene=gen, EG=j)) 
  }
  pos <- pos + 1
  if(pos %% 1000 == 0){
    print(pos)
  }
}
Sys.time() - time

eg_list <- new_df$EG
ref <- refseq[as.character(eg_list)]

time <- Sys.time()
df <- data.frame(RefSeq = character())
pos <- 1
for(i in unname(ref)){
  gen <- as.character(new_df$Gene[pos])
  for(j in i){
    if(j %like% "NM"){
      df <- rbind.fill(df, data.frame(RefSeq=j, gbm_dataset[gbm_dataset$Gene == gen, 2:ncol(gbm_dataset)]))
    }
  }
  pos <- pos + 1
  if(pos %% 1000 == 0){
    print(pos)
  }
}
Sys.time() - time

saveRDS(df, file="./df-reduced.rds")
df_copy <- readRDS("./df-reduced.rds")

Reduce(intersect, list(make.names(gbm_treatment_drug$bcr_patient_barcode), as.character(colnames(df))))

#subtype_patients <- gbm.subtypes_df[gbm.subtypes_df$Original.Subtype=="Mesenchymal", ]
subtype_patients <- gbm.subtypes_df
subtype_data <- df[, colnames(df) %in% make.names(subtype_patients$patient)]
subtype_data <- cbind(df[, 1:6], subtype_data)
#subtype_data <- subtype_data[!duplicated(subtype_data[, 2:ncol(subtype_data)]), ]
rownames(subtype_data) <- c(1:nrow(subtype_data))

control_group_names <- c("con.1", "con.2", "con.3", "con.4", "con.5")
gbm_group_names <- c()

for(i in c(1:(ncol(subtype_data)-6))){
  gbm_group_names <- c(gbm_group_names, paste("gbm.", i, sep=""))
}

colnames(subtype_data) <- c("RefSeq", control_group_names, gbm_group_names)

subtype_data_ref_list <- subtype_data$RefSeq
subtype_data$RefSeq <- NULL
factor_ind <- sapply(subtype_data, is.factor)
subtype_data[factor_ind] <- lapply(subtype_data[factor_ind], function(x) as.numeric(as.character(x)))
subtype_data <- cbind(subtype_data_ref_list, subtype_data)
colnames(subtype_data)[1] <- "RefSeq"

tcga_expset <- new("ExpressionSet", expr=as.matrix(subtype_data[, -1]))
fData(tcga_expset)[["gene_id"]] <- subtype_data[,1]
pData(tcga_expset)[["sample_name"]] <- colnames(subtype_data)[-1]

tcga_mirage_results <- MiRaGE(tcga_expset, species="HS")
tcga_mirage_results$P0[order(tcga_mirage_results$P0[,2])[1:20], ]

write.csv(tcga_mirage_results$P0[order(tcga_mirage_results$P0[,2])[1:20], ], "./MiRaGE_results/gbm_up.csv", row.names=F)
write.csv(tcga_mirage_results$P1[order(tcga_mirage_results$P1[,2])[1:20], ], "./MiRaGE_results/gbm_down.csv", row.names=F)
