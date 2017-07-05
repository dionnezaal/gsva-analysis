import sys
import pandas as pd

# Function that transforms the excel expression file to a expression table readable by python and R
def create_expression_table(FPKM):
	excel = pd.read_excel(FPKM, index_col=0)
	excel.columns = excel.columns.str.split('.').str[0] # Remove PD1 after patient

	# Create dictonairy of samples vs patients for columnname replacement
	sample_vs_patient = pd.read_table("sample_vs_patient_id.txt", index_col=1, sep=" ")
	dict_svsp = sample_vs_patient.to_dict()
	# Replace patient identifiers with the sample identifier in columns
	only_sample = dict_svsp.get('Accession')
	excel.rename(columns = only_sample, inplace=True)
	excel.index.name = 'Hugo_Symbol'
	return(excel)

def main():
	FPKM_file = sys.argv[1]

	# Transform excel expression file to expression table readable by python and R
	# and contains columns with sample identifiers instead of patient identifiers
	expression = create_expression_table(FPKM_file)
	# Write to file to be processed in R
	expression.to_csv("data_expression_file.txt", sep="\t")

main()