#!/bin/bash

###############################################################################
### Script to transform files from GEO (GSE78220) and supplemental material 
### from the melanoma Hugo (2016) paper (PMID:26997480) to import in cBioPortal

## Data files used:
## - GEO: Expression file: GSE78220_PatientFPKM.xlsx
## - GEO: Series Matrix File: GSE78220_series_matrix.txt.gz
## - GEO: SRA information table: sra_result.csv
## - Supplemental material: Table S1: mmc1.xls
###############################################################################

# paths to files
current_path=$(pwd)
expression_data="original_files/GSE78220_PatientFPKM.xlsx"
# Unzip the series matrix file with: gunzip GSE78220_series_matrix.txt.gz
series_matrix="original_files/GSE78220_series_matrix.txt"
sra_table="original_files/sra_result.csv"
supplemental_table1="original_files/mmc1.xls"

# Create new folder in current directory to put new files in
mkdir mela_hugo_2016

## First create study meta file
echo "type_of_cancer: mel
cancer_study_identifier: mela_hugo_2016
name: Melanoma (Hugo et al., 2016)
description: mRNA expressions in pre-treatment melanomas undergoing anti-PD-1 checkpoint inhibition therapy
citation: Genomic and Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic Melanoma, Hugo et al. Cell 2016
pmid: 26997480
short_name: Mela (Hugo)
add_global_case_list: true" > mela_hugo_2016/meta_study.txt

## Add (meta) cancer type files(for standalone docker)
echo "genetic_alteration_type: CANCER_TYPE
datatype: CANCER_TYPE
data_filename: cancer_type.txt" > mela_hugo_2016/meta_cancer_type.txt

echo -e "mel\tMelanoma\tmelanoma,melanoma metastatic\tBlack\ttissue" > mela_hugo_2016/cancer_type.txt

# Make meta patients and sample clinical data file
echo "cancer_study_identifier: mela_hugo_2016
genetic_alteration_type: CLINICAL
datatype: PATIENT_ATTRIBUTES
data_filename: data_clinical_patients.txt" > mela_hugo_2016/meta_clinical_patients.txt

echo "cancer_study_identifier: mela_hugo_2016
genetic_alteration_type: CLINICAL
datatype: SAMPLE_ATTRIBUTES
data_filename: data_clinical_samples.txt" > mela_hugo_2016/meta_clinical_samples.txt

## Create sample and patient clinical data files
# First create a stacked file with clinical attributes for each sample 
# identifier

# Grab all samples identifiers from series matrix
cat $series_matrix | grep "sample_id" | tr -d '"' | awk '{$1=""; print $0}' \
| tr " " "\n"  | sed '/^$/d' > sample_IDs.txt
# Grab all possible clinical attributes (sample characteristics in series 
# matrix file)
cat $series_matrix | grep "Sample_characteristics" | tr "\t" "\n" \
| sed '/Sample_characteristics_ch1/d' > sample_characteristics.txt

# Duplicate sample identifiers to be able to create stacked file
num_charac=$(cat $series_matrix | grep "characteristics" | wc -l)

# Make sure that there will be written to an empty file
echo -n > column_sample_id.txt

# Write sample identifiers multiple times to file
for ((i=1; i<=$num_charac; i++))
do 
	cat sample_IDs.txt >> column_sample_id.txt
done

# create stacked file by pasting the sample identifiers and clinical 
# attributes next to each other
paste -d "\t" column_sample_id.txt sample_characteristics.txt | tr -d "," > stacked_table.txt

# Remove empty clinical attributes for some sample identifiers
python scripts/drop_empty_attributes.py

# Sort stacked file on sample identifiers
cat nonan.csv | tr ":" "\t" | sort -k 2 > sorted_stacked_list.txt

# Create the clinical attribute files
python scripts/create_clinical_attribute_files.py 

## Create expression files
# meta expression file log2
echo "cancer_study_identifier: mela_hugo_2016
genetic_alteration_type: MRNA_EXPRESSION
datatype: CONTINUOUS
stable_id: rna_seq_mrna
show_profile_in_analysis_tab: false
profile_name: Log2 FPKM values
profile_description: Log2 FPKM values RNA-seq data
data_filename: log2_expression_file.txt" > mela_hugo_2016/meta_log2_expression_file.txt

# meta zscore expression file 
echo "cancer_study_identifier: mela_hugo_2016
genetic_alteration_type: MRNA_EXPRESSION
datatype: Z-SCORE
stable_id: rna_seq_mrna_median_Zscores
show_profile_in_analysis_tab: false
profile_name: Z-scores
profile_description: Z-scores RNA-seq data
data_filename: zscore_expression_file.txt" > mela_hugo_2016/meta_zscore_expression_file.txt

# Match sample and patient to transform in expression data
cat $sra_table | tr "," "\t" | awk '{print $2, $3}' | tr -d '":;' > sample_vs_patient_id.txt

# Create expression data table
python scripts/create_expression_data_file.py $expression_data

# Transform the expression table to log2 and zscore values
./scripts/get_log2_and_zscore_expression.R

## Create mutation information file from supplemental table
echo "cancer_study_identifier: mela_hugo_2016
genetic_alteration_type: MUTATION_EXTENDED
datatype: MAF
stable_id: mutations
show_profile_in_analysis_tab: true
profile_description: Mutation data
profile_name: Mutations
data_filename: mutations.txt" > mela_hugo_2016/meta_mutations.txt

python scripts/create_mutations_file.py $supplemental_table1

## Remove some files made during process
rm sample_IDs.txt sample_characteristics.txt sorted_stacked_list.txt sample_vs_patient_id.txt data_expression_file.txt nonan.csv stacked_table.txt column_sample_id.txt

# Print command to user that files are created (will also print this if there are errors)
echo "All files for import in cBioPortal are created and can be found in folder: mela_hugo_2016"