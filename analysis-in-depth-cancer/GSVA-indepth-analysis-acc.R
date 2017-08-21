##########################################################
## In depth cancer analysis of Adrenocortical carcinoma ##
## Author: Dionne Zaal                                  ##
##########################################################

# Package to draw volcano and boxplot
library(ggplot2)
# Package to reshape data frame
library(reshape2)

# Read in GSVA scores and patient data 
gsva_scores_acc <- read.table("/home/dionne/GSVA/Results_acc_tcga/GSVA_acc_tcga/acc_tcga_gsva_scores.txt", sep="\t", header=TRUE, row.names=1)
patient_info <- read.table("/home/dionne/GSVA/Results_acc_tcga/acc/tcga/data_bcr_clinical_data_patient.txt", sep="\t", header=TRUE)

# Make sure columnnames are the same as patient names in clinical file for matching between data
colnames(gsva_scores_acc) <- gsub("\\.", "-", colnames(gsva_scores_acc))
colnames(gsva_scores_acc) <- substr(colnames(gsva_scores_acc),1,nchar(colnames(gsva_scores_acc))-3)

# Grab only the patients which have GSVA score from the clinical file
patients_with_data <- patient_info[which(patient_info$PATIENT_ID %in% colnames(gsva_scores_acc)),]

# Sort the GSVA score columns in the same order as the clinical file
gsva_scores_acc <- gsva_scores_acc[,as.character(patients_with_data$PATIENT_ID)]  

# Grab the patients with and without tumor separatly
pt_with_tumor <- gsva_scores_acc[,patients_with_data$TUMOR_STATUS == "WITH TUMOR"]
pt_without_tumor <- gsva_scores_acc[,patients_with_data$TUMOR_STATUS == "TUMOR FREE"]

## Make volcano plot to show differences
# Will display the mean difference in GSVA score on the x axis and the -log10 pvalue of the difference between groups
mean_tumor <- apply(pt_with_tumor, 1, mean)
mean_without_tumor <- apply(pt_without_tumor, 1, mean)
diff_mean <- mean_tumor - mean_without_tumor

# Perform wilcox test to test if difference between groups is significant
wilcox_test_results <- apply(gsva_scores_acc, 1, function(x) wilcox.test(x[patients_with_data$TUMOR_STATUS == "WITH TUMOR"], x[patients_with_data$TUMOR_STATUS == "TUMOR FREE"])$p.value)

# Adjust the p-values with the holm-bonferroni method
adj_p_values_holm <- p.adjust(wilcox_test_results, method="holm")  # 500 significant differences

## Make volcano plot with all results
dat <- data.frame(difference_mean=diff_mean, p_value=-log10(wilcox_test_results))
# Color specific points in volcano plot because these are discussed in the results
G2_M <- c("REACTOME_G2_M_CHECKPOINTS", "REACTOME_G2_M_DNA_DAMAGE_CHECKPOINT", "SA_G2_AND_M_PHASES")
steroids <- c("GO_RESPONSE_TO_CORTICOSTEROID", "REACTOME_STEROID_HORMONES", "KEGG_STEROID_HORMONE_BIOSYNTHESIS")
p53 <- c("HALLMARK_P53_PATHWAY", "KEGG_P53_SIGNALING_PATHWAY", "BIOCARTA_P53_PATHWAY")
west <- c("WEST_ADRENOCORTICAL_TUMOR_UP", "WEST_ADRENOCORTICAL_CARCINOMA_VS_ADENOMA_UP", "WEST_ADRENOCORTICAL_TUMOR_DN", "WEST_ADRENOCORTICAL_CARCINOMA_VS_ADENOMA_DN")
west_marker <- c("WEST_ADRENOCORTICAL_TUMOR_MARKERS_UP", "WEST_ADRENOCORTICAL_TUMOR_MARKERS_DN")
sign_results <- c("REACTOME_IL_7_SIGNALING","REACTOME_SYNTHESIS_OF_BILE_ACIDS_AND_BILE_SALTS_VIA_24_HYDROXYCHOLESTEROL","KUMAMOTO_RESPONSE_TO_NUTLIN_3A_DN","KEGG_PRIMARY_BILE_ACID_BIOSYNTHESIS","REACTOME_BETA_DEFENSINS", "REACTOME_DEFENSINS", "GO_BILE_ACID_BIOSYNTHETIC_PROCESS","REACTOME_BILE_ACID_AND_BILE_SALT_METABOLISM", "REACTOME_SYNTHESIS_OF_BILE_ACIDS_AND_BILE_SALTS", "GO_DNA_POLYMERASE_COMPLEX", "INAMURA_LUNG_CANCER_SCC_UP", "MYLLYKANGAS_AMPLIFICATION_HOT_SPOT_1", "BIOCARTA_RANMS_PATHWAY","chr10q")
all <- c(G2_M, steroids, p53, west, west_marker,sign_results)
# Give names to the points in the plot to show in legend and match to
color_variable <- ifelse(adj_p_values_holm > 0.05, "Not significant", "Significant after correction")  
color_variable[G2_M] <- "G2/M phase"
color_variable[steroids] <- "Steroids"
color_variable[p53] <- "p53"
color_variable[west] <- "West"
color_variable[west_marker] <- "West marker"
color_variable[sign_results] <- "Highlighted significant results"

# Draw volcano plot
ggplot()+
  geom_point(data=dat[which(!rownames(dat)%in%all & adj_p_values_holm >0.05),], aes(x=difference_mean,y=p_value, colour=color_variable[which(!rownames(dat)%in%all & adj_p_values_holm >0.05)]), alpha=0.5) + 
  geom_point(data=dat[which(!rownames(dat)%in%all & adj_p_values_holm <= 0.05),], aes(x=difference_mean,y=p_value, colour=color_variable[which(!rownames(dat)%in%all & adj_p_values_holm <=0.05)]), alpha=0.5) + 
  geom_point(data=dat[rownames(dat)%in%G2_M,], aes(x=difference_mean,y=p_value, colour=color_variable[rownames(dat)%in%G2_M]), size=2) + 
  geom_point(data=dat[rownames(dat)%in%steroids,], aes(x=difference_mean,y=p_value, colour=color_variable[rownames(dat)%in%steroids]), size=2) + 
  geom_point(data=dat[rownames(dat)%in%p53,], aes(x=difference_mean,y=p_value, colour=color_variable[rownames(dat)%in%p53]), size=2) + 
  geom_point(data=dat[rownames(dat)%in%west,], aes(x=difference_mean,y=p_value, colour=color_variable[rownames(dat)%in%west]), size=2) + 
  geom_point(data=dat[rownames(dat)%in%west_marker,], aes(x=difference_mean,y=p_value, colour=color_variable[rownames(dat)%in%west_marker]), size=2) + 
  geom_point(data=dat[rownames(dat)%in%sign_results,], aes(x=difference_mean,y=p_value, colour=color_variable[rownames(dat)%in%sign_results])) + 
  geom_hline(yintercept=1.3, color="red") +
  scale_colour_manual(values=c("Significant after correction"="black", "Not significant"="darkgrey", "G2/M phase"="#0000FFFF", "Steroids"="#FFFF00FF", "p53" = "#00FF00FF", "West"="#00FFFFFF", "West marker"="#FFDB00FF", "Highlighted significant results"="#FF00FFFF"), breaks=c("Significant after correction", "Not significant", "G2/M phase", "Steroids", "p53", "West", "West marker",  "Highlighted significant results")) +
  labs(title="GSVA results - TCGA adrenocortical carcinoma",x="Difference in mean GSVA score", y="-log10 pvalue",  col="")+
  theme(plot.title=element_text(hjust=0.5)) +
  xlim(low=-0.6, high=0.6)

# Create boxplots, change p53 to show boxplot from other groups of gene sets
results <- cbind(pt_with_tumor[p53,], pt_without_tumor[p53,])
cancer_type <- rep(c(rep("with tumor", ncol(pt_with_tumor)), rep("without tumor", 
                                                                 ncol(pt_without_tumor))), length(p53))

# Results must be melted to be shown properly in the boxplot
melt_results <- melt(t(results))
melt_results$cancer_type <- cancer_type
colnames(melt_results) <- c("sample_id", "geneset", "GSVAscore", "neoplasm_status")

# Draw boxplot with within en between group differences can be seen
within_group <- ggplot(aes(y=GSVAscore, x=geneset, fill=neoplasm_status), data=melt_results) + geom_boxplot()
within_group + scale_fill_manual(values = c('lightgreen', 'darkgreen')) + labs(x=NULL, y="GSVA score", title="GSVA score distributions for p53 gene sets and neoplasm status") + theme(legend.title=element_blank(), plot.title=element_text(hjust=0.5)) +
  theme(legend.text=element_text(size=10)) + theme(legend.key.size =  unit(0.5, "in")) + ylim(low=-0.7, high=0.7)
