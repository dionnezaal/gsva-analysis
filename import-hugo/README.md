# Content mela_hugo_2016
Data from Melanoma patients from the research from Hugo et al. (2016) published in Cell (PMID:26997480). 

## folder mela_hugo_2016
Contains data files for import in cBioPortal. Study contains:
- Clinical data for patients and samples
- Expression data: log2 FPKM and Zscores
- Mutation data

## folder original files
Contains original files downloaded from the Gene Expression Omnibus or as supplemental information from the paper. These original files were used to create the data files in mela_hugo_2016. 

## folder scripts
Contains scripts to create the data files in mela_hugo_2016 from the original_files. 

## report.html
HTML file with validation report after importing the hugo script on local docker instance of cBioPortal

## transform_files_mela_hugo_2016.sh
This script is the main script that calls all scripts in the scripts folder to create the data files for cBioPortal from the original files. 

# Create data files cBioPortal from original files
When you want to create the cBioPortal data files again from the original analysis you should run: 

`./transform_files_mela_hugo_2016.sh`

This should be run in the same folder as transform_files_mela_hugo_2016.sh is located. The location of the original files are set in the top of the main script (directed to  original_files folder) and can be changed if necessary. 

