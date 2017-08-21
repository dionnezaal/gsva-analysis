#############################################################################
## Pan-cancer analysis of 31 TCGA datasets available on cBioPortal datahub ##
## Author: Dionne Zaal                                                     ##
#############################################################################

## Load packages for analysis
library("qusage") # Read gene set file
# Packages to create network of cancer studies
library(ggplot2)
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
tissues <- c(rep("gray",ncol(repr_scores)))
color_tissues <- c("gray", "#FFA500", "#DCDCDC", "#800080", "#808080")
tissues[c(which(colnames(repr_scores) %in% c("kich_tcga", "kirc_tcga", "kirp_tcga")))] <- "#FFA500"
tissues[c(which(colnames(repr_scores) %in% c("luad_tcga", "lusc_tcga")))] <- "#DCDCDC"
tissues[c(which(colnames(repr_scores) %in% c("pcpg_tcga", "acc_tcga")))] <- "#800080"
tissues[c(which(colnames(repr_scores) %in% c("gbm_tcga", "lgg_tcga")))] <- "#808080"

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

# Plot heatmap all genesets in the C2 collection
heatmap.2(as.matrix(repr_scores), labRow="", labCol=strsplit(colnames(repr_scores), "_tcga"), margins=c(8,3), cexCol=1.5, cexRow=1, col=hmcols, trace="none", dendrogram="none", Rowv=TRUE, Colv=TRUE)

# Plot heatmap selected genesets
heatmap.2(as.matrix(repr_scores_UPDOWN[UPandDOWN,]), labCol=strsplit(colnames(repr_scores), "_tcga"), cexCol=1.5, margins=c(8,35), cexRow=1, col=hmcols, trace="none", dendrogram="row", RowSideColors=types_UPanDOWN)

# Plot legends and paste later in heatmap
plot.new()
legend("topright",
       legend=c("Immune responses", "Cell Cycle", "Antigrowth signals", "Separate studies", "Metastasis", "Mutation"),
       col=c(color_types[1:6]),
       title=expression(bold("Gene set groups")),
       title.adj=0.20,
       pch=15,
       cex=1.2,
       pt.cex=1.5)

### Check if there are also tissue specific gene sets
# check with venn diagrams how many overlap there is between the same tissue types, with heatmap 
# with all studies we can see which genesets might be tissue specific

## kidney

kidney_up <- UP_genesets_per_cancertype[c("kich_tcga","kirc_tcga","kirp_tcga")]
kidney_down <- DOWN_genesets_per_cancertype[c("kich_tcga","kirc_tcga","kirp_tcga")]
venn_kidney_up <- venn(kidney_up)
overlap_up_kidney <- attributes(venn_kidney_up)$intersections$`kich_tcga:kirc_tcga:kirp_tcga`
venn_kidney_down <- venn(kidney_down)
overlap_down_kidney <- attributes(venn_kidney_down)$intersections$`kich_tcga:kirc_tcga:kirp_tcga`

# Remove some gene sets for more clear heatmap
overlap_up_kidney <- overlap_up_kidney[which(overlap_up_kidney %in% c("BIOCARTA_ABSCELL_PATHWAY", "BIOCARTA_BLYMPHOCYTE_PATHWAY", "REACTOME_BETA_DEFENSINS", "REACTOME_DEFENSINS", "BIOCARTA_IL5_PATHWAY", "BUDHU_LIVER_CANCER_METASTASIS_UP", "GALIE_TUMOR_STEMNESS_GENES", "SCHUHMACHER_MYC_TARGETS_DN", "REACTOME_ENDOSOMAL_VACUOLAR_PATHWAY","CHASSOT_SKIN_WOUND", "BRUNEAU_SEPTATION_ATRIAL","ALONSO_METASTASIS_EMT_DN","NAKAMURA_ALVEOLAR_EPITHELIUM","REACTOME_DIGESTION_OF_DIETARY_CARBOHYDRATE","SEIKE_LUNG_CANCER_POOR_SURVIVAL"))]

# heatmap all
tissues <- c(rep("gray",ncol(repr_scores)))
tissues[c(which(colnames(repr_scores) %in% c("kich_tcga", "kirc_tcga", "kirp_tcga")))] <- "#FFA500"
heatmap.2(as.matrix(repr_scores[c(overlap_up_kidney, overlap_down_kidney),sort(colnames(repr_scores))]), labCol=strsplit(colnames(repr_scores), "_tcga"), cexCol=1.5, margins=c(8,35), cexRow=1, col=hmcols, trace="none", ColSideColors=tissues)

## brain
brain_up <- UP_genesets_per_cancertype[c("gbm_tcga", "lgg_tcga")]
brain_down <- DOWN_genesets_per_cancertype[c(c("gbm_tcga", "lgg_tcga"))]
venn_brain_up <- venn(brain_up)
overlap_up_brain <- attributes(venn_brain_up)$intersections$`gbm_tcga:lgg_tcga`
venn_brain_down <- venn(brain_down)
overlap_down_brain <- attributes(venn_brain_down)$intersections$`gbm_tcga:lgg_tcga`

tissues <- c(rep("gray",ncol(repr_scores)))
tissues[c(which(colnames(repr_scores) %in% c("gbm_tcga", "lgg_tcga")))] <- "#808080"
# Heatmap with all overlapping gene sets
heatmap.2(as.matrix(repr_scores[c(overlap_up_brain, overlap_down_brain),sort(colnames(repr_scores))]), cexCol=1.5, margins=c(8,35), cexRow=1, col=hmcols, trace="none", ColSideColors=tissues)

# Smaller heatmap to emphasize the results
heatmap.2(as.matrix(repr_scores[c("REACTOME_ACETYLCHOLINE_NEUROTRANSMITTER_RELEASE_CYCLE","REACTOME_NOREPINEPHRINE_NEUROTRANSMITTER_RELEASE_CYCLE","REACTOME_DOPAMINE_NEUROTRANSMITTER_RELEASE_CYCLE","REACTOME_GABA_SYNTHESIS_RELEASE_REUPTAKE_AND_DEGRADATION"),sort(colnames(repr_scores))]), labCol=strsplit(colnames(repr_scores), "_tcga"), cexCol=1.5, margins=c(8,35), cexRow=1, col=hmcols, trace="none", ColSideColors=tissues)

## lung
lung_up <- UP_genesets_per_cancertype[c("luad_tcga", "lusc_tcga")]
lung_down <- DOWN_genesets_per_cancertype[c("luad_tcga", "lusc_tcga")]
venn_lung_up <- venn(lung_up)
overlap_up_lung <- attributes(venn_lung_up)$intersections$`luad_tcga:lusc_tcga`
venn_lung_down <- venn(lung_down)
overlap_down_lung <- attributes(venn_lung_down)$intersections$`luad_tcga:lusc_tcga`

tissues <- c(rep("gray67",ncol(repr_scores)))
tissues[c(which(colnames(repr_scores) %in% c("luad_tcga", "lusc_tcga")))] <- color_tissues[3]
heatmap.2(as.matrix(repr_scores[c(overlap_up_lung, overlap_down_lung),sort(colnames(repr_scores))]), cexCol=1.5, margins=c(8,35), cexRow=1, col=hmcols, trace="none", ColSideColors=tissues)

## Check if we can find significant differences between the kidney samples and all other samples when looking at all genesets
wilcox_test_results <- apply(repr_scores, 1, function(x) wilcox.test(x[c("kich_tcga", "kirc_tcga", "kirp_tcga")], x[!(colnames(repr_scores) %in% c("kich_tcga", "kirc_tcga", "kirp_tcga"))])$p.value)
adj_p_values_holm <- p.adjust(wilcox_test_results, method="holm")  #0 significant results
length(which(adj_p_values_holm < 0.05))

mean_kidney <- apply(repr_scores[,c("kich_tcga", "kirc_tcga", "kirp_tcga")], 1, mean)
mean_no_kidney <- apply(repr_scores[,!(colnames(repr_scores) %in% c("kich_tcga", "kirc_tcga", "kirp_tcga"))], 1, mean)
diff_mean <- mean_kidney - mean_no_kidney
dat <- data.frame(difference_mean=diff_mean, p_value=-log10(wilcox_test_results))

ggplot()+
  geom_point(data=dat[which(wilcox_test_results >0.05),], aes(x=difference_mean,y=p_value), colour="gray45", alpha=0.5) + 
  geom_point(data=dat[which(wilcox_test_results <= 0.05),], aes(x=difference_mean,y=p_value), colour="black", alpha=0.5)

## brain
wilcox_test_results <- apply(repr_scores, 1, function(x) wilcox.test(x[c("lgg_tcga", "gbm_tcga")], x[!(colnames(repr_scores) %in% c("lgg_tcga", "gbm_tcga"))])$p.value)
adj_p_values_holm <- p.adjust(wilcox_test_results, method="holm")  # 0 significant differences
length(which(adj_p_values_holm < 0.05))

mean_brain <- apply(repr_scores[,c("lgg_tcga", "gbm_tcga")], 1, mean)
mean_no_brain <- apply(repr_scores[,!(colnames(repr_scores) %in% c("lgg_tcga", "gbm_tcga"))], 1, mean)
diff_mean <- mean_brain - mean_no_brain
dat <- data.frame(difference_mean=diff_mean, p_value=-log10(wilcox_test_results))

ggplot()+
  geom_point(data=dat[which(wilcox_test_results >0.05),], aes(x=difference_mean,y=p_value), colour="gray45", alpha=0.5) + 
  geom_point(data=dat[which(wilcox_test_results <= 0.05),], aes(x=difference_mean,y=p_value), colour="black", alpha=0.5)


## lung
wilcox_test_results <- apply(repr_scores, 1, function(x) wilcox.test(x[c("luad_tcga", "lusc_tcga")], x[!(colnames(repr_scores) %in% c("luad_tcga", "lusc_tcga"))])$p.value)
adj_p_values_holm <- p.adjust(wilcox_test_results, method="holm")  # 0 significant differences
length(which(adj_p_values_holm < 0.05))

mean_lung <- apply(repr_scores[,c("luad_tcga", "lusc_tcga")], 1, mean)
mean_no_lung <- apply(repr_scores[,!(colnames(repr_scores) %in% c("luad_tcga", "lusc_tcga"))], 1, mean)
diff_mean <- mean_lung - mean_no_lung
dat <- data.frame(difference_mean=diff_mean, p_value=-log10(wilcox_test_results))

ggplot()+
  geom_point(data=dat[which(wilcox_test_results >0.05),], aes(x=difference_mean,y=p_value), colour="gray45", alpha=0.5) + 
  geom_point(data=dat[which(wilcox_test_results <= 0.05),], aes(x=difference_mean,y=p_value), colour="black", alpha=0.5)


## adrenal gland
wilcox_test_results <- apply(repr_scores, 1, function(x) wilcox.test(x[c("pcpg_tcga", "acc_tcga")], x[!(colnames(repr_scores) %in% c("pcpg_tcga", "acc_tcga"))])$p.value)
adj_p_values_holm <- p.adjust(wilcox_test_results, method="holm")  # 0 significant differences
length(which(adj_p_values_holm < 0.05))

mean_ag<- apply(repr_scores[,c("pcpg_tcga", "acc_tcga")], 1, mean)
mean_no_ag <- apply(repr_scores[,!(colnames(repr_scores) %in% c("pcpg_tcga", "acc_tcga"))], 1, mean)
diff_mean <- mean_ag - mean_no_ag
dat <- data.frame(difference_mean=diff_mean, p_value=-log10(wilcox_test_results))

ggplot()+
  geom_point(data=dat[which(wilcox_test_results >0.05),], aes(x=difference_mean,y=p_value), colour="gray45", alpha=0.5) + 
  geom_point(data=dat[which(wilcox_test_results <= 0.05),], aes(x=difference_mean,y=p_value), colour="black", alpha=0.5)

