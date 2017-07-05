import pandas as pd


# Function that removes empty attribute lines
def remove_empty_attributes(file):
	# Empy attributes ("") are recognized as NA in pandas
	series_matrix =pd.read_table(file, sep="\t", header=None)
	# Drop row if there ia a NA in this row
	nonan = series_matrix.dropna()
	return(nonan)

# Main function
def main():
	# Create stacked table without empty attribute lines
	nonan = remove_empty_attributes("stacked_table.txt")

	# Write to csv file
	nonan.to_csv("nonan.csv", sep="\t")

main()