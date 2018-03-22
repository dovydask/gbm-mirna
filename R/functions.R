read_exp_data <- function(pheno_file, exp_file, class_column){
  gbm_class_data <- read.delim(pheno_file)
  gbm_expr_data <- read.delim(exp_file)
  gbm_expr_data["Subtype"] <- gbm_class_data[class_column]
  gbm_expr_data <- gbm_expr_data[order(gbm_expr_data$Subtype),]
  rownames(gbm_expr_data) <- gbm_expr_data$Sample
  gbm_expr_data$Sample <- NULL
  return(gbm_expr_data)
}

linear_regression <- function(gbm_expr_data, max_pvalue = 0.025){
  data_results <- data.table(Gene=character(), pValue=numeric())
  i <- 1
  class_levels <- c(1, 2, 3, 4)
  classes <- as.numeric(gbm_expr_data$Subtype)
  for (gene in colnames(gbm_expr_data)) {
    print(i)
    i <- i + 1
    if (gene != "Subtype") {
      for (permutation in permn(class_levels)) {
        gene_expressions <- gbm_expr_data[[gene]]
        gen_df <- data.frame(classes, gene_expressions)
        gen_df$classes <- as.numeric(factor(gen_df$classes, levels=permutation))
        gen_df <- gen_df[order(gen_df$classes), ]
        
        model <- lm(gen_df$gene_expressions ~ gen_df$classes)
        pvalue <- as.numeric(anova(model)$'Pr(>F)'[1])
        if (pvalue < max_pvalue) {
          data_results <- rbindlist(list(data_results, list(gene, pvalue)))
        }
      }
    }
  }
  return(data_results)
}

filter_data <- function(gbm_expr_data, df, cutoff = 0.001){
  filtered <- df
  filtered <- filtered[order(filtered$gen_vec, filtered$pval_vec),]
  filtered <- subset(filtered, !duplicated(gen_vec))
  filtered <- filtered[filtered$pval_vec < 0.001, ]
  gbm_expr_data_filtered <- subset(gbm_expr_data, select=c(as.vector(filtered$gen_vec), "Subtype"))
  return(gbm_expr_data_filtered)
}

gene_plot <- function(gbm_exp_data, gene){
  classes <- as.numeric(gbm_expr_data$Subtype)
  gene_expressions <- gbm_expr_data[[gene]]
  plot(classes, gene_expressions, main="Gene Scatterplot")
  linear_model <- lm(gene_expressions ~ classes)
  abline(linear_model)
}