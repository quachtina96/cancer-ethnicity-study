import csv
import pandas as pd
import numpy as np
import pickle

data = np.load(open('../data/BRCA/BRCA.txt.matrix.npy'))
patients = []

with open('../data/BRCA/BRCA.txt.genes.csv', 'rb') as csvfile:
	reader = csv.reader(csvfile)
	count = 0
	for row in reader:
		for item in row:
			count += 1

with open('../data/BRCA/BRCA.txt.patients.csv', 'rb') as f:
	reader = csv.reader(f)
	patient_count = 0
	for row in reader:
		for item in row:
			patient_count += 1
			patients.append('-'.join(item.strip().split('-')[:3]))

assert patient_count == len(data[0])
assert count == data.shape[0]

clinical = pickle.load(open('../data/BRCA/BRCA.clinical/BRCA.clin.merged.picked.txt.saved.p'))

expression_data = list(data[0])
patient_data = {
		'expression': None, # use index of current gene
		'age': [],
		'gender': [],
		'race': [],
		'ethnicity': [],
		'stage': []
}

for i, patient in enumerate(patients):
	if patient in clinical:
		info = clinical[patient]
		patient_data['age'].append(info['years_to_birth'])
		patient_data['gender'].append(info['gender'])
		patient_data['race'].append(info['race'])
		patient_data['ethnicity'].append(info['ethnicity'])
		patient_data['stage'].append(info['pathologic_stage'])
	else:
		patients.remove(patient)
		del expression_data[i]

patient_data['expression'] = expression_data
df = pd.DataFrame(patient_data, index=patients)

# Now you have your data frame for gene 0!!!!!
