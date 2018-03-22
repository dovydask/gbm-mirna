# Setup -------------------------------------------------------------------

library(data.table)
library(ggplot2)
library(combinat)
library(furniture)
library(devtools)
library(ggbiplot)
library(rpart)
library(pracma)
library(factoextra)
library(FactoMineR)
library(party)
library(relaimpo)
library(earth)
library(Boruta)

setwd("C:/Users/Dovydas/Desktop/gbm/R/gene_analysis/Gene Analysis/TCGA/")
source("../functions.R")
pheno <- "../data/2018-02-09_TCGA_GBM_pheno.txt"
exp <- "../data/HG-U133A_2018-02-09_TCGA_GBM_expression.txt"
list_length <- 1000

data <- read_exp_data(pheno, exp, "Subtype_Verhaak_2010")
data <- data[complete.cases(data),]
#write.csv(data, file="./data.csv")

data <- readRDS("../../../mirna_analysis/TCGAWorkflow/expression_data.rds")
gbm_subtypes <- readRDS("../../../mirna_analysis/TCGAWorkflow/gbm_subtypes_df.rds")

data <- data[order(data$new_names), ]
data <- data[data$new_names != "?", ]
data <- data[!duplicated(data$new_names), ]
new_names <- unname(as.character(data$new_names))
patients <- colnames(data)
rownames(data) <- new_names
data$new_names <- NULL

data <- as.data.frame(t(data))

factors <- sapply(data, is.factor)
data[factors] <- lapply(data[factors], function(x) as.numeric(as.character(x)))
typeof(data[1, ncol(data)])

data <- cbind(data, gbm_subtypes[gbm_subtypes$patient %in% rownames(data), 2])
colnames(data)[ncol(data)] <- "Subtype"
data <- data[complete.cases(data), ]
data <- data[, -7570]

# Linear Regression -------------------------------------------------------

data.lr <- linear_regression(data)
saveRDS(data.lr, file="./data_lr.rds")
data.lr <- readRDS(file="./data_lr.rds")

# PCA ---------------------------------------------------------------------

pvalue_cutoff <- 0.000025
train.prop <- 0.9
pc_count <- 50

colnames(data.lr) <- c("Gene", "pValue")
data.ordered <- data.lr[order(data.lr$Gene, data.lr$pValue), ]
data.unique <- subset(data.ordered, !duplicated(Gene))
data.filtered <- data.unique[data.unique$pValue < pvalue_cutoff, ]
data.filtered <- data.filtered[order(data.filtered$pValue), ]
data.subtyped <- subset(data, select=c(as.vector(data.filtered$Gene), "Subtype"))
colnames(data.subtyped) <- make.names(colnames(data.subtyped))
#data.subtyped <- subset(data, select=c(as.vector(data.filtered$Gene), "HPRT1", "TBP", "GAPDH", "ACTB", "YWHAZ", "Subtype"))
data.complete <- data.subtyped[complete.cases(data.subtyped), ]
saveRDS(data.complete, file="./TCGA_data_complete.rds")
data.randomized <- data.complete[sample(nrow(data.complete)), ]
write(data.filtered$Gene[1:list_length], file="./lm_list.txt")

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

test <- data.classical[(train.prop * nrow(data.classical)):nrow(data.classical),]
test <- rbind(test, data.mesenchymal[(train.prop * nrow(data.mesenchymal)):nrow(data.mesenchymal),])
test <- rbind(test, data.neural[(train.prop * nrow(data.neural)):nrow(data.neural),])
test <- rbind(test, data.proneural[(train.prop * nrow(data.proneural)):nrow(data.proneural),])

library(randomForest)
model_rf <- randomForest(Subtype ~ ., data = train)
preds <- predict(model_rf, test)
table(preds)

sum((unname(preds) == test$Subtype) == TRUE)/length(test$Subtype)

pca.data <- train[, 1:ncol(train)-1]
subtypes <- train[, ncol(train)]
pca <- prcomp(pca.data, scale.=T, center=T)
pca.plot <- ggbiplot(pca, obs.scale = 1, var.scale = 1, 
              groups = subtypes, ellipse = TRUE, 
              circle = TRUE, var.axes = F)
pca.plot <- pca.plot + scale_color_discrete(name = '')
pca.plot <- pca.plot + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(pca.plot)

std_dev <- pca$sdev
var <- std_dev^2
prop_varex <- var/sum(var)
sum(prop_varex[1:pc_count])

train.data <- data.frame(Subtype = train$Subtype, pca$x)
train.data <- train.data[, 1:pc_count]
rpart.model <- rpart(Subtype ~ ., data = train.data, method = "anova")
test.data <- predict(pca, newdata = test)
test.data <- as.data.frame(test.data)
test.data <- test.data[,1:pc_count]
rpart.prediction <- predict(rpart.model, test.data)
rpart.prediction <- data.frame(as.list(round(rpart.prediction)))
rpart.prediction <- t(rpart.prediction)
rpart.prediction <- data.frame(rpart.prediction)
colnames(rpart.prediction) <- c("Subtype")
test.results <- test[match(rownames(rpart.prediction), make.names(rownames(test))), "Subtype"]
test.results <- as.numeric(test.results)
length(which(test.results == rpart.prediction$Subtype))/length(test.results)

# Most Important Variables ------------------------------------------------

loadings <- pca$rotation[, 1]
row_one_data <- scale(pca.data)[1, ]
calculated_pca <- sum(loadings * row_one_data)
all.equal(calculated_pca, pca$x[1, 1])

res.pca <- PCA(pca.data, graph = FALSE)
eigenvalues <- res.pca$eig
head(eigenvalues[, 1:2])
fviz_screeplot(res.pca, ncp=10)
contribution <- fviz_contrib(res.pca, choice = "var", axes = 1:2)
contribution.data <- contribution$data
contribution.data <- contribution.data[order(-contribution.data$contrib), ]
contribution.data

contribution.list <- sort(contribution.data$name[1:list_length])
write(as.character(contribution.list), file="./contr_list.txt")

cf1 <- cforest(Subtype ~ ., data=data.randomized, control=cforest_unbiased(mtry=2, ntree=50))
cf1
var.imp <- varimp(cf1)
var.imp <- sort(var.imp)
var.imp.auc <- varimpAUC(cf1)
var.imp.auc <- sort(var.imp.auc)
Reduce(intersect, list(names(var.imp[var.imp > 0]), names(var.imp.auc[var.imp.auc > 0])))
varimp.list <- var.imp[(length(var.imp)-list_length):length(var.imp)]
write(sort(names(varimp.list)), file="./varimp_list.txt")

boruta_output_numeric <- Boruta(as.numeric(Subtype) ~ ., data=data.randomized, doTrace=2)
boruta_output_numeric$finalDecision
boruta.list <- as.data.frame(boruta_output_numeric$finalDecision)
setDT(boruta.list, keep.rownames=TRUE)[]
colnames(boruta.list) <- c("Gene", "Status")
boruta.list.rejected <- boruta.list[boruta_output_numeric$finalDecision == "Rejected"]
boruta.list <- boruta.list[!(boruta.list$Gene %in% boruta.list.rejected$Gene)]
boruta.list <- boruta.list[boruta.list$Status == "Confirmed"]
write(boruta.list$Gene, file="./boruta_list.txt")

# Variance Plots ----------------------------------------------------------

plot(prop_varex, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")

plot(cumsum(prop_varex), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b")

# Gene Lists --------------------------------------------------------------

library(org.Hs.eg.db)

setwd("C:/Users/Dovydas/Desktop/gbm/R/gene_analysis/Gene Analysis")
sandmann_verhaak.list <- scan("./gene_lists/sandmann_verhaak.txt", what=character())
verhaak.list <- scan("./gene_lists/verhaak.txt", what=character())
lab.list <- scan("./gene_lists/neurooncolab.txt", what=character())

agilent.boruta <- scan("./Agilent-4502A/boruta_list.txt", what=character())
agilent.varimp <- scan("./Agilent-4502A/varimp_list.txt", what=character())
agilent.lm <- scan("./Agilent-4502A/lm_list.txt", what=character())

hg.boruta <- scan("./HG-U133A/boruta_list.txt", what=character())
hg.varimp <- scan("./HG-U133A/varimp_list.txt", what=character())
hg.lm <- scan("./HG-U133A/lm_list.txt", what=character())

tcga.boruta <- scan("./TCGA/boruta_list.txt", what=character())
tcga.varimp <- scan("./TCGA/varimp_list.txt", what=character())
tcga.lm <- scan("./TCGA/lm_list.txt", what=character())

hg.overlap <- Reduce(intersect, list(hg.lm, hg.varimp))
agilent.overlap <- Reduce(intersect, list(agilent.lm, agilent.varimp))
tcga.overlap <- Reduce(intersect, list(tcga.lm, tcga.varimp))
both.overlap <- c(hg.overlap, agilent.overlap)

write.csv(tcga.overlap, "./tcga-overlap.txt")
write.csv(all.overlap, "./all-overlap.txt")

all.overlap <- c(hg.overlap, agilent.overlap, tcga.overlap)
all.overlap

final_list <- sort(unique(all.overlap))
write.csv(final_list, "./final_list.txt", row.names=F, col.names=F)
genes.id <- select(org.Hs.eg.db, final_list, c("ENSEMBL", "ENTREZID"), "ALIAS")

write.csv(genes.id, file="./genes.csv", row.names=FALSE)

gene_overlap_list <- scan("../gene_overlap_list.txt", what=character())

data.overlap_list <- subset(data.complete, select=c(as.vector(gene_overlap_list), "Subtype"))

c(Reduce(intersect, list(hg.overlap, sandmann_verhaak.list)), Reduce(intersect, list(agilent.overlap, sandmann_verhaak.list)))
