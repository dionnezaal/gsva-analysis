#!/usr/bin/Rscript
## Perform GSVA analysis
start_time <- Sys.time()
cat("\n\n---> Load R libraries to perform the analysis: parallel, snow, qusage and GSVA\n\n")
# Load necessary libraries
library("parallel")  # To be able to do parallel analysis
library("snow")  # To perform parallel analysis
library("qusage") # To read in geneset file (.gmt)
library("GSVA") # To perform GSVA analysis

# Get arguments from command: 
# [1] name (and path) of expression file
# [2] name (and path) of geneset file
# [3] number of cores to use for the analysis
# [4] prefix for name outputfile
c_args <- commandArgs(TRUE)

expr_file <- c_args[1]
meta_expression_file <- c_args[2]
geneset_file <- c_args[3]
n_cores <- as.numeric(c_args[4])
n_bootstrap <- as.numeric(c_args[5])
prefix_out <- c_args[6]

## Load in data
# Expression file is expected to be separated with tabs, have samplenames as header 
# and rownames (Entrez Gene Ids) are available in the first column
cat(paste0("\n\n---> Load expression file ", expr_file, " and geneset file ",  geneset_file, "\n\n"))
expr <- read.delim(expr_file, sep="\t", header=T, row.names=2, quote="")
expr <- expr[,-1]
expr_norm <- log2(expr + 1)
mexp <- rowMeans(expr_norm)
expr_norm_high <- expr_norm[mexp > 1, ]

# Geneset file should be in gmt format
genesets <- read.gmt(geneset_file)  # Will put genesets in named list

## Perform GSVA analysis
cat(paste0("---> Perform GSVA analysis in parallel, will use ", n_cores, " cores for calculation. Number of bootstraps set to ", n_bootstrap, "\n\n"))
gsva_result <- gsva(as.matrix(expr_norm_high), genesets, method="gsva", no.bootstraps=n_bootstrap, parallel.sz=n_cores, parallel.type="SOCK")
gsva_scores <- data.frame("geneset_id" = rownames(gsva_result$es.obs), gsva_result$es.obs, check.names = F)
gsva_bootstrap <- data.frame("geneset_id" = rownames(gsva_result$bootstrap$p.vals.sign),gsva_result$bootstrap$p.vals.sign, check.names = F)

colnames(gsva_scores) <- gsub("\\.", "-", colnames(gsva_scores))
colnames(gsva_bootstrap) <- gsub("\\.", "-", colnames(gsva_bootstrap))

## Write gsva scores to file
write.table(gsva_scores, paste0(prefix_out, "_gsva_scores.txt"), quote = F, sep = "\t", col.names = T, row.names = F)
write.table(gsva_bootstrap, paste0(prefix_out, "_gsva_pvalues.txt"), quote = F, sep = "\t", col.names = T, row.names = F)

cat(paste0("\n\n---> Output files written to ", prefix_out, "_gsva_scores.txt and ", prefix_out, "_gsva_pvalues.txt\n\n"))

cat(paste0("\n\n---> Create meta datafiles for gsva scores and pvalues, also create new case list"))

### Write meta files
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
geneset_def_version: msigdb_v6.0_and_cbio")
write(meta_scores, paste0(prefix_out, "_meta_gsva_scores.txt"))

meta_pvalues <- paste0("cancer_study_identifier: ", study_id, "
genetic_alteration_type: GENESET_SCORE
datatype: P-VALUE
stable_id: gsva_pvalues
source_stable_id: gsva_scores
profile_name: GSVA p-values
profile_description: P-values calculated with GSVA boostrapping method (n=", n_bootstrap, ").
data_filename: ", tail(strsplit(prefix_out, "/")[[1]],1), "_gsva_pvalues.txt
geneset_def_version: msigdb_v6.0_and_cbio")
write(meta_pvalues, paste0(prefix_out, "_meta_gsva_pvalues.txt"))

case_list <- paste0(c("cancer_study_identifier: ", study_id, "
stable_id: ", study_id, "_gsva_scores
case_list_name: Tumor Samples with GSVA data
case_list_description: All samples with GSVA data
case_list_category: all_cases_with_gsva_data
case_list_ids:    ", paste0(colnames(gsva_scores)[2:ncol(gsva_scores)], collapse = "\t")), collapse = "")
write(case_list, paste0(prefix_out, "_cases_GSVA.txt"))

cat(paste0("\n\n---> Meta files written to ", prefix_out, "_meta_gsva_scores.txt, ", prefix_out, "_meta_gsva_pvalues.txt and ", prefix_out, "_case_list.txt\n\n"))

save.image(file=paste0(prefix_out, "_GSVA_calc.RData"))

cat(paste0("\n\n---> R environment variables written to ", prefix_out, "_GSVA_calc.RData\n\n"))

end_time <- Sys.time()
cat(paste0("Elapsed time: ", difftime(end_time, start_time, units="mins"), " mins\n\n"))
