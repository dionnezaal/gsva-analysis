#!/usr/bin/env Rscript
# Read in expression data made with create_expression_data_file.py
hugo_expr <- read.table("data_expression_file.txt", header=TRUE)

# Normalizes expression by taking the log2 to minimize outliers
hugo_expr_norm <- log2(hugo_expr[,-1] + 1)

# Make table which writes the correct expression file for cBioPortal
# Genes need to be in separate column to be able to write index name Hugo_Symbol to file
Hugo_Symbol <- hugo_expr[,1]
rownames(hugo_expr_norm) <- NULL
writable_hugo_expr_norm <- cbind(Hugo_Symbol, hugo_expr_norm)
write.table(writable_hugo_expr_norm, file="mela_hugo_2016/log2_expression_file.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

# Calculate zscore: (expression value - mean gene) / standard deviation gene
zscore_matrix <- t(apply(hugo_expr_norm, 1, function(x) ((x-mean(x)) / sd(x))))

# Write zscore also to expression file
colnames(zscore_matrix) <- colnames(hugo_expr_norm)
rownames(zscore_matrix) <- NULL
writable_zscores_matrix <- cbind(Hugo_Symbol,as.data.frame(zscore_matrix))
write.table(writable_zscores_matrix, file="mela_hugo_2016/zscore_expression_file.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
