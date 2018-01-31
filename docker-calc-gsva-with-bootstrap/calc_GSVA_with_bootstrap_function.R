#!/usr/bin/Rscript

###############################################################
## Script to perform GSVA analysis for studies in datahub    ##
## Input: expression file (rows: genes, columns: samples),   ##
##		    meta expression file, geneset file (in gmt),       ##
##		    number of cores to use for analysis, number of     ##
## 		    bootstraps to use for analysis, prefix for output  ##
## 		    files (with desired path)                          ##
## Output: GSVA score and p-value file (genesets in rows,    ##
##		     samples in columns), meta score and p-value file, ##
##         case list file and the Rdata from the analysis	   ##
## Author: Dionne Zaal									                     ##
###############################################################

cat("\n\n---> Load R libraries to perform the analysis: parallel, snow, qusage and GSVA\n\n")
library("parallel")  # To be able to do parallel analysis
library("snow")  # To perform parallel analysis
library("qusage") # To read in geneset file (.gmt)
library("GSVA") # To perform GSVA analysis

# Get arguments from command: 
# [1] name (and path) of expression file
# [2] name (and path) of geneset file
# [3] name (and path) of meta expression file
# [4] prefix for name outputfile
# [5] number of cores to use for the analysis
# [6] number of bootstraps for the analysis
c_args <- commandArgs(TRUE)
expr_file <- c_args[1]
geneset_file <- c_args[2]
meta_expression_file <- c_args[3]
prefix_out <- c_args[4]
n_cores <- as.numeric(c_args[5])
n_bootstrap <- as.numeric(c_args[6])

# Load and normalize expression file
cat(paste0("\n\n---> Load expression file ", expr_file, " and geneset file ",  geneset_file, "\n\n"))
expr <- read.delim(expr_file, sep="\t", header=T, row.names=2, quote="")
expr <- expr[,-1]
expr_norm <- log2(expr + 1)
mexp <- rowMeans(expr_norm)
expr_norm_high <- expr_norm[mexp > 1, ]

# Load genesets, geneset file should be in gmt format
genesets <- read.gmt(geneset_file)  # Will put genesets in named list

# Calculate original gene set scores
cat("\n\n---> Calculate GSVA scores original genesets with ", n_cores," number of cores\n\n")
full_set_gsva_result <- gsva(as.matrix(expr_norm_high), genesets, method="gsva", parallel.sz=n_cores, parallel.type="SOCK")

# Write scores to file
cat(paste0("\n\n---> Output file original GSVA written to ", prefix_out, "_gsva_scores.txt \n\n"))
colnames(full_set_gsva_result$es.obs) <- gsub("\\.", "-", colnames(full_set_gsva_result$es.obs))
write.table(full_set_gsva_result$es.obs, paste0(prefix_out, "_gsva_scores.txt"), quote = F, sep = "\t", col.names = T, row.names = T)

# If the user indicated that bootstraps should be done, do bootstrapping
if (n_bootstrap > 0){
  cat("\n\n---> Randomize gene sets and calculate bootstrap scores")
  # Create empty list for new genesets
  new_genesets <- vector("list", length(genesets))
  names(new_genesets) <- names(genesets)
  
  # Create list to add all bootstrap scores separatly
  list_scores <- vector("list", n_bootstrap)
  names(list_scores) <- paste0("bootstrap_", 1:n_bootstrap)
  all_genes_geneset <- unique(as.vector(unlist(genesets)))
  
  # source function to resample gene sets
  source("func_sampling_same_dist.R")
  
  for (i in 1:n_bootstrap){
    # Resample gene sets from the same distribution with function in separate R file
    new_genesets <- sample_geneset_from_dist(genesets)
    
    # method that does not use the same distribution but grabs random sets:
    # all_genes_geneset <- unique(as.vector(unlist(genesets)))
    # for (geneset in 1:length(genesets)){
    # new_genesets[geneset][[1]] <- c(sample(all_genes_geneset, length(genesets[geneset][[1]]), replace=FALSE))
    # }
    
    # Calculate GSVA scores for bootstrapped gene sets
    gene_resamp_gsva <- gsva(as.matrix(expr_norm_high), new_genesets, method="gsva", parallel.sz=n_cores, parallel.type="SOCK")
    results <- data.frame(matrix(NA, nrow=nrow(full_set_gsva_result$es.obs), ncol=ncol(full_set_gsva_result$es.obs)))
    
    # Save 'bootstrap' scores to list with dataframes (size of dataframes are the same and therefore can complete be added to list)
    rownames(results) = rownames(full_set_gsva_result$es.obs)
    colnames(results) = colnames(full_set_gsva_result$es.obs)
    results[rownames(gene_resamp_gsva$es.obs),] <- gene_resamp_gsva$es.obs
    list_scores[[i]] <- results
  }

  # Create empty dataframe to add p-values
  pvalues_calc <- data.frame(matrix(NA, nrow=nrow(full_set_gsva_result$es.obs), ncol=ncol(full_set_gsva_result$es.obs)))
  rownames(pvalues_calc) = rownames(full_set_gsva_result$es.obs)
  colnames(pvalues_calc) = colnames(full_set_gsva_result$es.obs)
  
  # Calculating p-values with bootstrapped values
  cat("\n\n---> Calculate pvalues")
  for (sample in 1:ncol(pvalues_calc)){
    for (geneset in 1:nrow(pvalues_calc)){
      # Get all scores from each sample and gene sets for the amount of bootstraps
      bootstrap_scores <- sapply(list_scores, '[', geneset, sample) 
      
      # Positive and negative scores needed separatly for calculation
      pos_scores <- bootstrap_scores[which(bootstrap_scores >= 0)]
      neg_scores <- bootstrap_scores[which(bootstrap_scores < 0)]
      if (full_set_gsva_result$es.obs[geneset, sample] >= 0) {
        pvalues_calc[geneset,sample] <- sum(pos_scores >= full_set_gsva_result$es.obs[geneset, sample]) / (length(pos_scores) +1)
      } else {
        pvalues_calc[geneset,sample] <- sum(neg_scores <= full_set_gsva_result$es.obs[geneset, sample]) / (length(neg_scores) +1)
      }

      # In case no scores are above (or below in case of negative) the original, set the pvalue to the lowest possible
      if (pvalues_calc[geneset,sample] == 0){
        pvalues_calc[geneset,sample] = 1/n_bootstrap
      }
    }
  }
  
  # Write p-value results to file
  cat(paste0("\n\n---> Output file bootstrap written to ", prefix_out, "_gsva_pvalues.txt\n\n"))
  colnames(pvalues_calc) <- gsub("\\.", "-", colnames(pvalues_calc))
  write.table(pvalues_calc, paste0(prefix_out, "_gsva_pvalues.txt"), quote = F, sep = "\t", col.names = T, row.names = F)
  
}


### Write meta files
cat(paste0("\n\n---> Create meta datafiles for gsva scores and pvalues, also create new case list"))
# Get versions R and GSVA
r_version <- as.character(getRversion())
gsva_version <- as.character(packageVersion("GSVA"))

# Read in expression meta file
meta_expression <- read.table(meta_expression_file, sep=":", header=FALSE)
study_id <- trimws(as.character(meta_expression[which(meta_expression[,1] == "cancer_study_identifier"), 2]))
source_stable_id <- trimws(as.character(meta_expression[which(meta_expression[,1] == "stable_id"), 2]))

meta_scores <- paste0("cancer_study_identifier: ", study_id, "
genetic_alteration_type: GENESET_SCORE
datatype: GSVA-SCORE
stable_id: gsva_scores
source_stable_id: ", source_stable_id, "
profile_name: GSVA scores
profile_description: GSVA scores for MSigDB v6.0 genesets calculated with GSVA version ", gsva_version,", R version ", r_version, "
data_filename: ", tail(strsplit(prefix_out, "/")[[1]],1), "_gsva_scores.txt
geneset_def_version: msigdb_v6.0")
write(meta_scores, paste0(prefix_out, "_meta_gsva_scores.txt"))

# In case bootstrapping is done make meta file which indicates the amount of
# bootstraps that were done, otherwise create a dummy pvalue file and
# indicate this in the meta p-value file
if (n_bootstrap > 0){
  meta_pvalues <- paste0("cancer_study_identifier: ", study_id, "
genetic_alteration_type: GENESET_SCORE
datatype: P-VALUE
stable_id: gsva_pvalues
source_stable_id: gsva_scores
profile_name: GSVA p-values
profile_description: P-values calculated with GSVA boostrapping method (n=", n_bootstrap, ").
data_filename: ", tail(strsplit(prefix_out, "/")[[1]],1), "_gsva_pvalues.txt
geneset_def_version: msigdb_v6.0")
  write(meta_pvalues, paste0(prefix_out, "_meta_gsva_pvalues.txt"))
} else {
  # create dummy p-values instead of empty file
  dummy_pvalues <- (full_set_gsva_result$es.obs * 0) + 0.01
  gsva_pvalues <- data.frame("geneset_id" = rownames(full_set_gsva_result$es.obs), dummy_pvalues, check.names = F)
  colnames(gsva_pvalues) <- gsub("\\.", "-", colnames(gsva_pvalues))
  write.table(gsva_pvalues, paste0(prefix_out, "_gsva_pvalues.txt"), quote = F, sep = "\t", col.names = T, row.names = F)
  
  # create also meta p-value file
  meta_pvalues <- paste0("cancer_study_identifier: ", study_id, "
genetic_alteration_type: GENESET_SCORE
datatype: P-VALUE
stable_id: gsva_pvalues
source_stable_id: gsva_scores
profile_name: GSVA p-values
profile_description: Dummy P-values, no bootstrap done.
data_filename: ", tail(strsplit(prefix_out, "/")[[1]],1), "_gsva_pvalues.txt
geneset_def_version: msigdb_v6.0")
  write(meta_pvalues, paste0(prefix_out, "_meta_gsva_pvalues.txt"))
}

case_list <- paste0(c("cancer_study_identifier: ", study_id, "
stable_id: ", study_id, "_gsva_scores
case_list_name: Tumor Samples with GSVA data
case_list_description: All samples with GSVA data
case_list_category: all_cases_with_gsva_data
case_list_ids:    ", paste0(colnames(full_set_gsva_result$es.obs), collapse = "\t")), collapse = "")
write(case_list, paste0(prefix_out, "_cases_GSVA.txt"))

cat(paste0("\n\n---> Meta files written to ", prefix_out, "_meta_gsva_scores.txt, ", prefix_out, "_meta_gsva_pvalues.txt and ", prefix_out, "_case_list.txt\n\n"))

save.image(file=paste0(prefix_out, "_GSVA_calc.RData"))

cat(paste0("\n\n---> R environment variables written to ", prefix_out, "_GSVA_calc.RData\n\n"))

