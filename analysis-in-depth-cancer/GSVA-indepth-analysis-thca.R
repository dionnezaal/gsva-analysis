################################################
## In depth cancer analysis of Thyroid cancer ##
## Author: Dionne Zaal                        ##
################################################

# Load packages
library("ggplot2") # For drawing box and volcanoplot
library("reshape2") # To melt the data for boxplot

# Load GSVA scores, sample and patient data
gsva_scores_thca <- read.table("/home/dionne/GSVA/Results_thca_tcga/import_thca_tcga/data_gsva_scores.txt", sep="\t", header=TRUE, row.names=1)
sample_info <- read.table("/home/dionne/GSVA/Results_thca_tcga/thca/tcga/data_bcr_clinical_data_sample.txt", sep="\t", header=TRUE)
patient_info <- read.table("/home/dionne/GSVA/Results_thca_tcga/thca/tcga/data_bcr_clinical_data_patient.txt", sep="\t", header=TRUE)

# first filter sample file with only samples we have
colnames(gsva_scores_thca) <- gsub("\\.", "-", colnames(gsva_scores_thca))
sub_sample_info <- sample_info[which(sample_info$SAMPLE_ID %in% colnames(gsva_scores_thca)),]
order_gsva_scores_thca <- gsva_scores_thca[,as.character(sub_sample_info$SAMPLE_ID)]

# Seperate tables of different cancer types
PTC_patients <- order_gsva_scores_thca[,which(sub_sample_info$CANCER_TYPE_DETAILED == "Papillary Thyroid Cancer")]
FTC_patients <- order_gsva_scores_thca[,which(sub_sample_info$CANCER_TYPE_DETAILED == "Follicular Thyroid Cancer")]

## Make volcano plot to show differences
# Will display the mean difference in GSVA score on the x axis and the -log10 pvalue of the difference between groups
mean_papillary <- apply(PTC_patients, 1, mean)
mean_follicular <- apply(FTC_patients, 1, mean)
diff_mean <- mean_papillary - mean_follicular

# Perform wilcox test to determine if difference is significant
wilcox_test_results <- apply(gsva_scores_thca, 1, function(x) wilcox.test(x[colnames(PTC_patients)], x[colnames(FTC_patients)])$p.value)

# Adjust p-value for multiple testing
adj_p_values_holm <- p.adjust(wilcox_test_results, method="holm")  # 7808 significant differences

# Data for volcanoplot
dat <- data.frame(difference_mean=diff_mean, p_value=-log10(wilcox_test_results))

# Define certain points in the volcano plot for coloring
nfkappab <- c("GILMORE_CORE_NFKB_PATHWAY", "BIOCARTA_NFKB_PATHWAY")
wnt <- c("BIOCARTA_WNT_PATHWAY", "REACTOME_SIGNALING_BY_WNT", "KEGG_WNT_SIGNALING_PATHWAY")
translocation <- c("chr2q13", "chr3p25")
other_studies <- c("DELYS_THYROID_CANCER_UP","FONTAINE_PAPILLARY_THYROID_CARCINOMA_UP","FONTAINE_FOLLICULAR_THYROID_ADENOMA_DN", "DELYS_THYROID_CANCER_DN", "FONTAINE_PAPILLARY_THYROID_CARCINOMA_DN", "FONTAINE_FOLLICULAR_THYROID_ADENOMA_UP")
sign_results <- c("GO_MITOCHONDRIAL_ELECTRON_TRANSPORT_CYTOCHROME_C_TO_OXYGEN","chr4p11", "ST_TYPE_I_INTERFERON_PATHWAY", "REACTOME_REGULATION_OF_COMPLEMENT_CASCADE", "HINATA_NFKB_TARGETS_KERATINOCYTE_DN", "GO_CELLULAR_GLUCURONIDATION", "GO_ODORANT_BINDING")
all <- c(other_studies,nfkappab,wnt,sign_results)

color_variable <- ifelse(adj_p_values_holm > 0.05, "Not significant", "Significant")  
color_variable[other_studies] <- "Other thyroid cancer studies"
color_variable[nfkappab] <- "Nf kappaB"
color_variable[wnt] <- "Wnt signaling"
color_variable[sign_results] <- "Highlighted significant results"

# Draw volcanoplot
ggplot()+
  geom_point(data=dat[which(!rownames(dat)%in%all & adj_p_values_holm >0.05),], aes(x=difference_mean,y=p_value, colour=color_variable[which(!rownames(dat)%in%all & adj_p_values_holm >0.05)]), alpha=0.5) + 
  geom_point(data=dat[which(!rownames(dat)%in%all & adj_p_values_holm <= 0.05),], aes(x=difference_mean,y=p_value, colour=color_variable[which(!rownames(dat)%in%all & adj_p_values_holm <=0.05)]), alpha=0.5) + 
  geom_point(data=dat[rownames(dat)%in%other_studies,], aes(x=difference_mean,y=p_value, colour=color_variable[rownames(dat)%in%other_studies]), alpha=0.8, size=2) + 
  geom_point(data=dat[rownames(dat)%in%nfkappab,], aes(x=difference_mean,y=p_value, colour=color_variable[rownames(dat)%in%nfkappab]), alpha=0.8, size=2) + 
  geom_point(data=dat[rownames(dat)%in%wnt,], aes(x=difference_mean,y=p_value, colour=color_variable[rownames(dat)%in%wnt]), alpha=0.8, size=2) + 
  geom_point(data=dat[rownames(dat)%in%sign_results,], aes(x=difference_mean,y=p_value, colour=color_variable[rownames(dat)%in%sign_results]), alpha=0.8) + 
  geom_hline(yintercept=1.3, color="red") +
  scale_colour_manual(values=c("Significant"="darkgrey", "Not significant"="black", "Other thyroid cancer studies"=rainbow(6)[6], "Nf kappaB"=rainbow(6)[3], "Wnt signaling"=rainbow(6)[2], "Translocation"=rainbow(6)[4], "Highlighted significant results"=rainbow(6)[5]), breaks=c("Significant", "Not significant", "Other thyroid cancer studies", "Nf kappaB", "Wnt signaling", "Translocation", "Highlighted significant results")) +
  labs(title="GSVA results - TCGA Thyroid cancer",x="Difference in mean GSVA score", y="-log10 pvalue",  col="")+
  theme(plot.title=element_text(hjust=0.5)) +
  xlim(low=-0.6, high=0.6)

## Draw boxplots of certain genesets
# replace other studies for different genesets
results <- cbind(pt_with_tumor[other_studies,], pt_without_tumor[other_studies,])
cancer_type <- rep(c(rep("PTC", ncol(PTC_patients)), rep("FTC", ncol(FTC_patients))), 
                   length(other_studies))

melt_results <- melt(t(results))
melt_results$cancer_type <- cancer_type
colnames(melt_translocation_results) <- c("sample_id", "geneset", "GSVAscore", "cancer_type")

within_group <- ggplot(aes(y=GSVAscore, x=geneset, fill=cancer_type), data=melt_translocation_results) + geom_boxplot()
within_group + scale_fill_manual(values = c('lightgreen', 'darkgreen')) + labs(x=NULL, y="GSVA score", title="GSVA score per geneset and cancer subtype") + theme(legend.title=element_blank(), plot.title=element_text(hjust=0.5))

## Check if more patients with papillary had radiotherapy after surgery
patient_info <- read.table("/home/dionne/GSVA/Results_thca_tcga/thca/tcga/data_bcr_clinical_data_patient.txt", sep="\t", header=TRUE)
sub_patient_info <- patient_info[which(patient_info$PATIENT_ID %in% sub_sample_info$PATIENT_ID),]

treatment_yes <- sub_patient_info[which(sub_patient_info$RADIATION_TREATMENT_ADJUVANT == "YES"),]
treatment_no <- sub_patient_info[which(sub_patient_info$RADIATION_TREATMENT_ADJUVANT == "NO"),]
treatment_NA <- sub_patient_info[which(sub_patient_info$RADIATION_TREATMENT_ADJUVANT == "[Not Available]"),] # No unavailable treatments

sum(treatment_yes$PATIENT_ID %in% strsplit(colnames(PTC_patients), "-01")) #134
sum(treatment_yes$PATIENT_ID %in% strsplit(colnames(FTC_patients), "-01")) #22
sum(treatment_no$PATIENT_ID %in% strsplit(colnames(PTC_patients), "-01")) # 58
sum(treatment_no$PATIENT_ID %in% strsplit(colnames(FTC_patients), "-01")) #18
sum(treatment_NA$PATIENT_ID %in% strsplit(colnames(PTC_patients), "-01")) # 202
sum(treatment_NA$PATIENT_ID %in% strsplit(colnames(FTC_patients), "-01")) # 66

perc_papillary <- (134/ (ncol(PTC_patients) - 202)) * 100  # 66.66%
perc_follicular <- (22/(ncol(FTC_patients) - 66)) * 100 # 55%

## Quality check for differential analysis
# Check if random groups also result in this in the same amount of differential genesets
number_sign_genes <- c()
number_non_sign_genes <- c()
for (i in 1:100){
  resample_names <- sample(colnames(gsva_scores_thca), ncol(gsva_scores_thca), replace=FALSE)
  group1_scores <- gsva_scores_thca[,resample_names[1:106]]
  group2_scores <- gsva_scores_thca[,resample_names[107:509]]
  
  mean_group1 <- apply(group1_scores, 1, mean)
  mean_group2 <- apply(group2_scores, 1, mean)
  diff_mean <- mean_group1 - mean_group2

  wilcox_test_results <- apply(gsva_scores_thca, 1, function(x) wilcox.test(x[colnames(group1_scores)], x[colnames(group2_scores)])$p.value)
  adj_pvalues <- p.adjust(wilcox_test_results, method="holm")
  number_sign_genes <- c(number_sign_genes, sum(adj_pvalues <= 0.05))
  number_non_sign_genes <- c(number_non_sign_genes, sum(adj_pvalues > 0.05))
}

# Cannot only plot from the last resampling
plot(diff_mean, -log10(wilcox_test_results), col=color_points, xlab="GSVA score difference", ylab="-log10 pvalue")

# But can determine max and minimum of significant genes
min(number_sign_genes)
max(number_sign_genes)
# These numbers are much less than original therefore the differential analysis done seems to be correct

# Density curve to check distributions of plots
melt_FTC <- melt(FTC_patients)
melt_PTC <- melt(PTC_patients)
ggplot() + geom_density(dat=melt_FTC, aes(x=value, color="FTC")) + 
  geom_density(dat=melt_PTC,aes(x=value, color="PTC")) + 
  labs(x="GSVA scores", y="Density", title="Distribution of GSVA scores - TCGA Thyroid Cancer") +
  scale_colour_manual(values=c("FTC" = "lightgreen", "PTC" = "darkgreen"), name="Density")

