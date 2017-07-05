#!/usr/bin/Rscript

###############################################################
## Calculate representative GSVA score per geneset per study ##
## Input: GSVA scores file (rows: genesets, columns: samples)##
##        prefix for name of output files (with path)        ##
## Output: file with representative GSVA per geneset         ##
## Author: Dionne Zaal                                       ##
###############################################################

c_args <- commandArgs(TRUE)
gsva_scores_file <- c_args[1]
prefix_output <- c_args[2]

gsva_scores <- read.table(gsva_scores_file, sep="\t", header=TRUE, row.names=1)

repr_geneset_scores <- c(rep(NA, nrow(gsva_scores)))
names(repr_geneset_scores) <- rownames(gsva_scores)

# Calculate a representative score for each geneset for this study
for (geneset in 1:nrow(gsva_scores)){
  ## Separate in positive and negative scores
  pos_scores <- gsva_scores[geneset, which(gsva_scores[geneset,] >= 0)]
  neg_scores <- gsva_scores[geneset, which(gsva_scores[geneset,] < 0)]
  
  # If there are no positive scores add the negative score quantile as representative score
  if (ncol(pos_scores) == 0){
    repr_geneset_scores[rownames(gsva_scores)[geneset]] <- quantile(as.matrix(neg_scores), probs=c(0.25))
  } 
  # If there is no negative score add the positive score quantile as representative score
  else if (ncol(neg_scores) == 0){
    repr_geneset_scores[rownames(gsva_scores)[geneset]] <- quantile(as.matrix(pos_scores), probs=c(0.75))
  } 
  # If there are both positive and negative calculate both scores and put in the highest as representative
  else {
    quantile_pos_scores <- quantile(as.matrix(pos_scores), probs=c(0.75))
    quantile_neg_scores <- quantile(as.matrix(neg_scores), probs=c(0.25))
    repr_geneset_scores[rownames(gsva_scores)[geneset]] <- ifelse(quantile_pos_scores > abs(quantile_neg_scores), quantile_pos_scores, quantile_neg_scores)
  }
}

# Write results to file
write.table(repr_geneset_scores, paste0(prefix_output, "_repr_GSVA_scores.txt"), quote=FALSE, sep="\t", row.names=names(repr_geneset_scores), col.names=FALSE)