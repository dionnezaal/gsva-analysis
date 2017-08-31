#!/usr/bin/Rscript

################################################
## Pan-Cancer analysis test                   ##
## Input: Expression values of specific TCGA  ##
##        studies                             ##
## Output: Message if the same result         ##
##         is found as the original for the   ##
##         pan-cancer analysis                ##
## Author: Dionne Zaal                        ##
################################################

#install.packages(c('snow', 'DBI', 'RSQLite', 'RCurl', 'xtable', 'openssl', 'Rmpi', 'covr', 'rlecuyer'), repos='http://cran-mirror.cs.uu.nl/', dependencies=TRUE)
#install.packages(c('qusage', 'GSEABase', 'GSVA'), repos='http://bioconductor.org/packages/3.5/bioc')
library("parallel")  # To be able to do parallel analysis
library("snow")  # To perform parallel analysis
library("qusage") # To read in geneset file (.gmt)
library("GSVA") # To perform GSVA analysis

# Get arguments from command: 
# [1] path to datahub files
# [2] number of cores to use for GSVA
# [3] name (and path) of gene set file with complete MSigDB collection
# [4] name (and path) of gene set file with only C2 collection MSigDB
# [5] TRUE or FALSE inidicating that GSVA scores should be calculated or not
c_args <- commandArgs(TRUE)
path_to_datahub <- c_args[1]  # Path to folders with datahub studies (public folder)
n_cores <- c_args[2]  # numbers of cores to use for GSVA analysis
geneset_file <- c_args[3] # Gene set file for GSVA should be all collection from msigdb
c2_geneset_file <- c_args[4] # For the pan cancer analysis only the c2 collection is used
calc_GSVA <- c_args[5] # In case GSVA still needs to be done this parameter should be set to TRUE, otherwise set to FALSE

## Studies used in the original analysis
studies <- c("meso_tcga","kich_tcga","ov_tcga","thym","dlbc_tcga","gbm_tcga","laml_tcga","chol_tcga", 
             "acc_tcga","lgg_tcga","kirp_tcga","pcpg_tcga","kirc_tcga","thca_tcga","prad_tcga","tgct_tcga",
             "esca_tcga","stad_tcga","luad_tcga","cesc_tcga","ucs_tcga","paad_tcga","uvm_tcga","lusc_tcga",
             "hnsc_tcga","skcm_tcga","sarc_tcga","blca_tcga","brca_tcga","lihc_tcga","coadread_tcga")

# Original gene sets found in common in at least 20 cancer studies
original_UPandDOWN <- c("BOWIE_RESPONSE_TO_TAMOXIFEN","WATANABE_ULCERATIVE_COLITIS_WITH_CANCER_DN",
                          "KEGG_ALLOGRAFT_REJECTION","FINETTI_BREAST_CANCER_KINOME_GREEN",
                          "ROETH_TERT_TARGETS_DN","BIOCARTA_DC_PATHWAY","BIOCARTA_CTLA4_PATHWAY",
                          "REACTOME_PD1_SIGNALING","REACTOME_PHOSPHORYLATION_OF_CD3_AND_TCR_ZETA_CHAINS",
                          "BIOCARTA_ASBCELL_PATHWAY","BIOCARTA_BLYMPHOCYTE_PATHWAY",
                          "BIOCARTA_NO2IL12_PATHWAY","REACTOME_TRANSLOCATION_OF_ZAP_70_TO_IMMUNOLOGICAL_SYNAPSE",
                          "CHAN_INTERFERON_PRODUCING_DENDRITIC_CELL","BIOCARTA_CTL_PATHWAY",
                          "BIOCARTA_THELPER_PATHWAY","BIOCARTA_TCRA_PATHWAY","BIOCARTA_TCYTOTOXIC_PATHWAY",
                          "BIOCARTA_TCAPOPTOSIS_PATHWAY","CHASSOT_SKIN_WOUND","SCHUHMACHER_MYC_TARGETS_DN",
                          "BIOCARTA_IL5_PATHWAY","ZHANG_INTERFERON_RESPONSE","BUDHU_LIVER_CANCER_METASTASIS_UP",
                          "BOWIE_RESPONSE_TO_EXTRACELLULAR_MATRIX","BIERIE_INFLAMMATORY_RESPONSE_TGFB1",
                          "MOSERLE_IFNA_RESPONSE","FARMER_BREAST_CANCER_CLUSTER_4","REACTOME_UNWINDING_OF_DNA",
                          "ZHAN_MULTIPLE_MYELOMA_PR_UP","MONTERO_THYROID_CANCER_POOR_SURVIVAL_UP","LY_AGING_MIDDLE_DN",
                          "FINETTI_BREAST_CANCER_KINOME_RED","CROSBY_E2F4_TARGETS","ZERBINI_RESPONSE_TO_SULINDAC_DN",
                          "LIANG_SILENCED_BY_METHYLATION_DN","KUMAMOTO_RESPONSE_TO_NUTLIN_3A_DN")

# Read geneset file, geneset file should be in gmt format
genesets <- read.gmt(geneset_file)  # Will put genesets in named list
c2_genesets <- read.gmt(c2_geneset_file)

# Big matrix to paste all representative scores per study
repr_scores <- data.frame(study=rep(NA,length(names(c2_genesets))), row.names=names(c2_genesets))

# For each studies calculate the GSVA score and representative score
for (i in 1:length(studies)){
  
  if (calc_GSVA == "TRUE"){
    expr_file <- paste0(path_to_datahub, studies[i], "/", gsub("_tcga", "", studies[i]),"/tcga/data_RNA_Seq_v2_expression_median.txt")
    
    ## Read and normalize expression file
    cat(paste0("\n\n---> Load expression file ", expr_file))
    expr <- read.delim(expr_file, sep="\t", header=T, row.names=2, quote="")
    expr <- expr[,-1]
    expr_norm <- log2(expr + 1)
    mexp <- rowMeans(expr_norm)
    expr_norm_high <- expr_norm[mexp > 1, ]
    
    cat(paste0("\n\n---> Calculate GSVA scores study ",studies[i]),"\n\n")
    # Calculate original gene set scores
    gsva_result <- gsva(as.matrix(expr_norm_high), genesets, method="gsva", parallel.sz=n_cores, parallel.type="SOCK")
    gsva_scores <- gsva_result$es.obs
  } else {
    gsva_file <- paste0(path_to_datahub, studies[i], "/", gsub("_tcga", "", studies[i]),"/tcga/gsva_scores.txt")
    gsva_scores <- read.table(gsva_file, row.names = 1, header=TRUE)
  }
  
  repr_geneset_scores <- c(rep(NA, nrow(gsva_scores)))
  names(repr_geneset_scores) <- rownames(gsva_scores)
  
  cat(paste0("\n\n---> Calculate representative scores study ",studies[i]),"\n\n")
  for (geneset in 1:nrow(gsva_scores)){
    ## Separate in positive and negative scores
    pos_scores <- gsva_scores[geneset, which(gsva_scores[geneset,] >= 0)]
    neg_scores <- gsva_scores[geneset, which(gsva_scores[geneset,] < 0)]
    
    # If there are no positive scores add the negative score quantile as representative score
    if (length(pos_scores) == 0){
      repr_geneset_scores[rownames(gsva_scores)[geneset]] <- quantile(as.matrix(neg_scores), probs=c(0.25))
    } 
    # If there is no negative score add the positive score quantile as representative score
    else if (length(neg_scores) == 0){
      repr_geneset_scores[rownames(gsva_scores)[geneset]] <- quantile(as.matrix(pos_scores), probs=c(0.75))
    } 
    # If there are both positive and negative calculate both scores and put in the highest as representative
    else {
      quantile_pos_scores <- quantile(as.matrix(pos_scores), probs=c(0.75))
      quantile_neg_scores <- quantile(as.matrix(neg_scores), probs=c(0.25))
      repr_geneset_scores[rownames(gsva_scores)[geneset]] <- ifelse(quantile_pos_scores > abs(quantile_neg_scores), quantile_pos_scores, quantile_neg_scores)
    }
  }
  repr_scores[,i] <- repr_geneset_scores[rownames(repr_scores)]
}

colnames(repr_scores) <- studies

## Select from the representative scores the X number of highest and lowest scoring genesets
## per cancer study/type
UP_genesets_per_cancertype <- list()
DOWN_genesets_per_cancertype <- list()
num_genesets <- 100

for (cancer_study in 1:ncol(repr_scores)){
  # Save the X number highest and lowest genesets (aka UP and DOWN)
  UP <- tail(rownames(repr_scores[order(repr_scores[,cancer_study]),]), n=num_genesets)
  DOWN <- head(rownames(repr_scores[order(repr_scores[,cancer_study]),]), n=num_genesets)
  
  # Save genesets in separate list for highest and lowest genesets for a cancer study/type
  name_study <- paste0(colnames(repr_scores[cancer_study]))
  UP_genesets_per_cancertype[[name_study]] <- UP
  DOWN_genesets_per_cancertype[[name_study]] <- DOWN
}

## Heatmap of all genesets that were identified as UP or DOWN in any of the cancer studies 
genesets_UPandDOWN <- unique(c(unlist(UP_genesets_per_cancertype, use.names=F), unlist(DOWN_genesets_per_cancertype, use.names=F)))

# Save only representative scores which were defined as UP and DOWN in most studies
repr_scores_UPDOWN <- repr_scores[genesets_UPandDOWN,]

# Find genesets that are up or down in at least half of the cancer studies
UP_count <- data.frame(count=rep(0, length(genesets_UPandDOWN)), row.names=genesets_UPandDOWN)
DOWN_count <- data.frame(count=rep(0, length(genesets_UPandDOWN)), row.names=genesets_UPandDOWN)

# Count how often certain genesets are found in the cancer studies
for (i in 1:length(genesets_UPandDOWN)){
  UP_count[i,1] <- length(which(unlist(UP_genesets_per_cancertype) %in% genesets_UPandDOWN[i]))
  DOWN_count[i,1] <- length(which(unlist(DOWN_genesets_per_cancertype) %in% genesets_UPandDOWN[i]))
}

# Save only gene sets which are found in more than 20 gene sets as UP or DOWN gene sets
UPandDOWN <- c(rownames(UP_count)[which(UP_count > 20)], rownames(DOWN_count)[which(DOWN_count > 20)])

# check the amount of differences between the vectors
if (length(setdiff(original_UPandDOWN, UPandDOWN)) == 0){
  cat("\n\nNo differences between original and current gene sets from pan-cancer analysis\n\n")
} else {
  cat("\n\nDifferences found between original and current gene sets from pan-cancer analysis\n\n")
  cat(paste0("Original gene set not in current gene sets: ", setdiff(original_UPandDOWN, UPandDOWN), "\n"))
  cat(paste0("Current gene set not in original gene sets: ", setdiff(UPandDOWN,original_UPandDOWN), "\n"))
}