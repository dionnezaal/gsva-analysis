#############################################################################
## Pan-cancer analysis of 31 TCGA datasets available on cBioPortal datahub ##
## Author: Dionne Zaal                                                     ##
#############################################################################

## Load packages for analysis
library("qusage") # Read gene set file
# Packages to create network of cancer studies
library(network)
library(sna)
library(ggplot2)
library(intergraph)
library(GGally)
# Packages for heatmap and colors in heatmap
library("RColorBrewer")
library("gplots")

# Use representative scores per study to do the cross cancer analysis
repr_scores <- read.table("/home/dionne/GSVA/representative_scores/all_repr_scores.txt", sep="\t", row.names=1)

# Select only curated genesets to limit the results
c2_genesets <- read.gmt("/home/dionne/Downloads/c2.all.v6.0.entrez.gmt")
repr_scores <- repr_scores[names(c2_genesets),]

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

## Calculate number of overlapping genesets for each combination of cancer studies
## separatly for the UP and DOWN genesets and display the overlap in a network

comb_cancertypes <- combn(colnames(repr_scores),2) # Get all pairs of cancer studies/types
# Appending amount of genesets in common
UPandDOWN_network <- data.frame(cancer_type1=as.character(), cancer_type2=as.character(), common_genesets=as.integer(), stringsAsFactors = FALSE)
# Network variable, a 1 will be appended if there is an edge between two cancer studies/types
network <- data.frame(matrix(ncol=ncol(repr_scores), nrow=ncol(repr_scores)))
rownames(network) <- colnames(repr_scores)
colnames(network) <- colnames(repr_scores)
# Set edges to own cancer study/type to zero
diag(network) <- 0

edge_sizes <- c()  # Set sizes for edge thickness between studies/types 

for (comb in 1:ncol(comb_cancertypes)){
  common_UP <- length(intersect(UP_genesets_per_cancertype[[comb_cancertypes[1,comb]]], UP_genesets_per_cancertype[[comb_cancertypes[2,comb]]]))
  common_DOWN <- length(intersect(DOWN_genesets_per_cancertype[[comb_cancertypes[1,comb]]], DOWN_genesets_per_cancertype[[comb_cancertypes[2,comb]]]))
  common_both <- as.integer(common_UP + common_DOWN)
  UPandDOWN_network[comb,] <- c(comb_cancertypes[1,comb], comb_cancertypes[2,comb], common_both)
  # Put edges in network
  # More than 20 UP and DOWN genesets should be in common between the gene sets before an edge is shown
  test_UP <- ifelse(common_UP > 20, 1, 0)
  test_DOWN <- ifelse(common_DOWN > 20, 1, 0)
  both <- test_UP + test_DOWN
  network[UPandDOWN_network$cancer_type1[comb], UPandDOWN_network$cancer_type2[comb]] <- ifelse(both == 2, 1, 0)
  # Determine edge thickness
  if (network[UPandDOWN_network$cancer_type1[comb], UPandDOWN_network$cancer_type2[comb]] == 1){
    edge_sizes <- c(edge_sizes, round((common_UP+common_DOWN)/50)+0.5)
  }
  
  # Give top right of the matrix a zero, so no edge will be drawn
  # Because the edge will already be drawn 
  network[UPandDOWN_network$cancer_type2[comb], UPandDOWN_network$cancer_type1[comb]] <- 0
}

net = network(network, directed=FALSE)

# Give tissue types to specific studies to later color nodes in network by tissue type
name_tissues <- c(rep("Distinct",ncol(repr_scores)))
name_tissues[c(which(colnames(repr_scores) %in% c("ov_tcga","ucs_tcga", "cesc_tcga")))] <- "Uterus"
name_tissues[c(which(colnames(repr_scores) %in% c("kich_tcga", "kirc_tcga", "kirp_tcga")))] <- "Kidney"
name_tissues[c(which(colnames(repr_scores) %in% c("luad_tcga", "lusc_tcga")))] <- "Lung"
name_tissues[c(which(colnames(repr_scores) %in% c("pcpg_tcga", "acc_tcga")))] <- "Adrenal gland"
name_tissues[c(which(colnames(repr_scores) %in% c("gbm_tcga", "lgg_tcga")))] <- "Brain"
name_tissues[c(which(colnames(repr_scores) %in% c("hnsc_tcga", "skcm_tcga")))] <- "Skin"

# Attach the tissue type to the network
net %v% "tissues" = name_tissues

# Colors for different tissuetypes
colors_palette <- c("gray67", colorRampPalette(c("mediumpurple2","cyan2"))(6))
names(colors_palette) <- c("Distinct", "Uterus", "Kidney", "Lung", "Adrenal gland", "Brain", "Skin")

# Plot network
ggnet2(net, label=gsub("_tcga", "", colnames(repr_scores)), color="tissues", palette=colors_palette, color.legend="Tissue type", label.size=3, edge.size=edge_sizes, edge.color = "black")

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

# Color columns in heatmap according to tissue type
tissues <- c(rep("gray67",ncol(repr_scores)))
color_tissues <- colorRampPalette(c("mediumpurple2","cyan2"))(6)
tissues[c(which(colnames(repr_scores) %in% c("ov_tcga","ucs_tcga", "cesc_tcga")))] <- color_tissues[1]
tissues[c(which(colnames(repr_scores) %in% c("kich_tcga", "kirc_tcga", "kirp_tcga")))] <- color_tissues[2]
tissues[c(which(colnames(repr_scores) %in% c("luad_tcga", "lusc_tcga")))] <- color_tissues[3]
tissues[c(which(colnames(repr_scores) %in% c("pcpg_tcga", "acc_tcga")))] <-color_tissues[4]
tissues[c(which(colnames(repr_scores) %in% c("gbm_tcga", "lgg_tcga")))] <- color_tissues[5]
tissues[c(which(colnames(repr_scores) %in% c("hnsc_tcga", "skcm_tcga")))] <- color_tissues[6]

# Define colors for GSVA scores in heatmap
hmcols<-colorRampPalette(c("limegreen","white","deeppink2"))(40)

# Give different colors to specific processes found with gene sets
color_types <- colorRampPalette(c("gold","firebrick3"))(6)
types_UPanDOWN <- c(rep(color_types[1], length(UPandDOWN)))
cell_cycle <- c("CHASSOT_SKIN_WOUND", "REACTOME_UNWINDING_OF_DNA","ZHAN_MULTIPLE_MYELOMA_PR_UP", "ZERBINI_RESPONSE_TO_SULINDAC_DN", "CROSBY_E2F4_TARGETS", "FINETTI_BREAST_CANCER_KINOME_RED")
types_UPanDOWN[which(UPandDOWN %in% cell_cycle)] <- color_types[2]
antigrowth <- c("SCHUHMACHER_MYC_TARGETS_DN")
types_UPanDOWN[which(UPandDOWN %in% antigrowth)] <- color_types[3]
sep_studies <- c("LY_AGING_MIDDLE_DN", "LIANG_SILENCED_BY_METHYLATION_DN","MONTERO_THYROID_CANCER_POOR_SURVIVAL_UP","WATANABE_ULCERATIVE_COLITIS_WITH_CANCER_DN")
types_UPanDOWN[which(UPandDOWN %in% sep_studies)] <- color_types[4]
metastasis <- c("BUDHU_LIVER_CANCER_METASTASIS_UP", "BOWIE_RESPONSE_TO_EXTRACELLULAR_MATRIX", "FARMER_BREAST_CANCER_CLUSTER_4")
types_UPanDOWN[which(UPandDOWN %in% metastasis)] <- color_types[5]
mutation <- c("KUMAMOTO_RESPONSE_TO_NUTLIN_3A_DN")
types_UPanDOWN[which(UPandDOWN %in% mutation)] <- color_types[6]

# Plot heatmap
heatmap.2(as.matrix(repr_scores_UPDOWN[UPandDOWN,]), labCol=strsplit(colnames(repr_scores), "_tcga"), cexCol=1.5, margins=c(8,35), cexRow=1, col=hmcols, trace="none", ColSideColors=tissues, dendrogram="both", RowSideColors=types_UPanDOWN)

# Plot legends and paste later in heatmap
plot.new()
legend("topright",
       legend=c("Immune responses", "Cell Cycle", "Antigrowth signals", "Separate studies", "Metastasis", "Mutation"),
       col=c(color_types[1:6]),
       title=expression(bold("Gene set groups")),
       title.adj=0.20,
       pch=15,
       cex=.8,
       pt.cex=1.5)

plot.new()
legend("topright",
       legend=c("Distinct", "Uterus", "Kidney", "Lung", "Adrenal gland", "Brain", "Skin"),
       col=c("gray67", color_tissues),
       title=expression(bold("Tissue types")),
       title.adj=0.20,
       pch=15,
       cex=.8,
       pt.cex=1.5)