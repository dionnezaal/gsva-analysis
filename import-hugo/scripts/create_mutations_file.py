import sys
import pandas as pd
import warnings

def create_mutation_table(suppl_table1):
	# Supplemental table 2 with sheet S1D contains mutation information
	supplemental = pd.read_excel(suppl_table1, sheetname="S1D", header=1)
	# Create new chromosome column which only has chromosome number or letter (X,Y)
	supplemental['Chromosome'] = supplemental['Chr'].str.strip('chr')
	# Select columns necessary for mutation table
	sub_supl = supplemental[['Gene','Sample','MutType', 'Aamut', 'Chromosome', 'Pos']]
	# Rename the columns
	sub_supl.columns = ['Hugo_Symbol', 'PatientID', 'Variant_Classification', 'HGVSp_Short', 'Chromosome', 'Start_Position']
	# In cBioPortal mutations are displayed per sample, therefore mutations from patient 27 needs to be duplicated
	sub_supl['Tumor_Sample_Barcode'] = sub_supl['PatientID'].replace(to_replace="Pt27", value="Pt27a") 
	mutation_data = sub_supl.append(sub_supl.loc[sub_supl['Tumor_Sample_Barcode'] == "Pt27a"].replace(to_replace="Pt27a", value="Pt27b"))
	
	# Now change the patient identiefiers to sample identifiers
	sample_data = pd.read_table("mela_hugo_2016/data_clinical_samples.txt", sep="\t", header=4)
	sample_data['PATIENT_ID'][19] = 'Pt27a'
	sample_data['PATIENT_ID'][20] = 'Pt27b'
	sample_data['PATIENT_ID'] = sample_data['PATIENT_ID'].str.lstrip() # Remove spaces for right match
	patient_vs_sample = sample_data.set_index('PATIENT_ID').to_dict()['SAMPLE_ID'] # Make dict of patient vs sample column
	# Add sample IDs as Tumor_Sample_Barcode to the mutation file
	mutation_data['Tumor_Sample_Barcode'] = mutation_data['Tumor_Sample_Barcode'].replace(patient_vs_sample)
	del mutation_data['PatientID']
	# From patients which do not have expression data also mutation information was available, remove these patients
	filter = mutation_data['Tumor_Sample_Barcode'].str.contains("GSM")
	mutation_data = mutation_data[filter]
	# Add 'p.' in front of mutation to be recognizes correctly
	mutation_data['HGVSp_Short'] = 'p.' + mutation_data['HGVSp_Short'].astype(str)
	mutation_data['HGVSp_Short'] = mutation_data['HGVSp_Short'].replace(to_replace="p.nan", value="NA")
	return(mutation_data)


def main():
	# Supress pandas warning SettingWithCopyWarning
	warnings.simplefilter(action='ignore', category=Warning)

	supplemental_table1 = sys.argv[1]
	# Transform supplemental table to mutation file for cBioPortal
	mutation_data = create_mutation_table(supplemental_table1)
	# Write mutation file to correct location
	mutation_data.to_csv("mela_hugo_2016/mutations.txt", sep="\t", index=False)

main()