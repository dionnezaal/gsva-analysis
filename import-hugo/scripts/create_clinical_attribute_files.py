import pandas as pd 

# Open stacked table and transform to table with samples as rows and attributes as columns
stacked = pd.read_table("sorted_stacked_list.txt", sep="\t")
stacked.columns = ['SampleID', 'variables', 'values']
table = stacked.pivot(index='SampleID', columns='variables', values='values')

# Retrieve different columns for sample and patient attribute file
sample_columns = ['patient id' ,'anatomical location', 'tissue', 'biopsy time']
sample_table = table[sample_columns]
patient_columns = ['patient id', 'age (yrs)', 'anti-pd-1 response', 'disease status', 'gender', 'overall survival (days)', 'previous mapki', 'study site','treatment', 'vital status']
patient_table = table[patient_columns]

# Delete duplicated patient lines (Pt27 has two samples and therefore a duplicated line)
patient_nodup = patient_table.drop_duplicates(subset='patient id')
# Set as index the patient identifier and remove the sample identifier as index
patient_nosamp = patient_nodup.reset_index()
patient_nosamp1 = patient_nosamp.drop('SampleID', 1)
patient_nosamp1.index = patient_nosamp1['patient id']
patient_nosamp1 = patient_nosamp1.drop('patient id', 1)

# Remove whitespaces in front and back of column values
patient_nosamp1[patient_nosamp1.columns] = patient_nosamp1.apply(lambda x: x.str.strip())

## Replace some values to more meaningful ones or ones cBioPortal acceptes
# Replace disease status with full name (found in paper: Dickson, Paxton V., and Jeffrey E. Gershenwald. "Staging and prognosis of cutaneous melanoma." Surgical oncology clinics of North America 20.1 (2011): 1-17)
values_replacement = {'M0': 'No distant metastases', 'M1a': 'Distant skin, subcutaneous, or nodal metastases', 'M1b': 'Lung metastases', 'M1c': 'All other visceral and distant metastasis'}
# Assumed that no distant metastased means disease free, other still have disease. No disease free days/months in data file
values_replacement_disease_free = {'M0': 'DiseaseFree', 'M1a': 'Recurred/Progressed', 'M1b': 'Recurred/Progressed', 'M1c': 'Recurred/Progressed'}
patient_nosamp1['disease free status'] = patient_nosamp1['disease status'].replace(values_replacement_disease_free)
patient_nosamp1['disease status'] = patient_nosamp1['disease status'].replace(values_replacement)
patient_nosamp1['previous mapki'] = patient_nosamp1['previous mapki'].replace({'Y':'Yes', 'N':'No'})
patient_nosamp1['vital status'] = patient_nosamp1['vital status'].replace({'Dead':'DECEASED', 'Alive':'LIVING'})

# Calculate survival in months instead of days
copy_column = patient_nosamp1['overall survival (days)'].replace({'NA': 0})
days = copy_column.tolist()
months = [int(i)/30 for i in days]
months = [x if x != 0 else 'NA' for x in months]  # Replace 0 with NA again
patient_nosamp1['overal survival (months)'] = months # add months column to table

# Add necessary lines for clinical attributes files
# Samle clinical attributes
top_lines = pd.DataFrame(index=['#Sample identification', '#STRING','#1', 'SAMPLE_ID'],columns=[sample_table.columns])
top_lines.index.name = '#SampleID'
top_lines.loc['#Sample identification'] = top_lines.columns # Extra info about columns, for now the same as the columns
top_lines.loc['#STRING'] = ['STRING', 'STRING', 'STRING', 'STRING']
top_lines.loc['#1'] = [1] * 4
top_lines.loc['SAMPLE_ID'] = ['PATIENT_ID','ANATOMICAL_LOCATION', 'TISSUE','TIME_BIOPSY']
# Now we combine the first lines (top_lines) with the table with full results (table)
sample_table.index.name = "#SampleID"
all_frames = [top_lines, sample_table]
result = pd.concat(all_frames)
# Write clinical sample table to file
result.to_csv("mela_hugo_2016/data_clinical_samples.txt", sep="\t")

# Patient clinical attributes
top_lines = pd.DataFrame(index=['#Patient identification', '#STRING','#1', 'PATIENT_ID'],columns=[patient_nosamp1.columns])
top_lines.index.name = '#PatientID'
top_lines.loc['#Patient identification'] = ['age in years', 'anti-pd-1 response', 'disease status', 'gender', 'overall survival in days', 'previous mapk inhibitor (resistance to treatment can occur)', 'study site','treatment', 'vital status', 'disease free status', 'overal survival days in months (calculated from days divide by 30)'] # Extra info about columns, for now the same as the columns
top_lines.loc['#STRING'] = ['NUMBER', 'STRING', 'STRING','STRING','NUMBER', 'STRING', 'STRING', 'STRING', 'STRING','STRING','NUMBER']
top_lines.loc['#1'] = [1] * 11
top_lines.loc['PATIENT_ID'] = ['AGE','RESPONSE_ANTIPD1', 'DS_STATUS', 'GENDER', 'SURVIVAL_DAYS', 'PREVIOUS_MAPKI', 'STUDY_SITE', 'TREATMENT', ' OS_STATUS', 'DFS_STATUS', 'OS_MONTHS']
# Now we combine the first lines (top_lines) with the table with full results (table)
patient_nosamp1.index.name = "#PatientID"
all_frames = [top_lines, patient_nosamp1]
result = pd.concat(all_frames)
# Write clinical patient table to file
result.to_csv("mela_hugo_2016/data_clinical_patients.txt", sep="\t")

