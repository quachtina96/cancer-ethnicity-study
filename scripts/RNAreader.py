import os
import csv
import numpy as np
import pickle

def read_rna(directory):
	print "Conducting Random Forest analysis on RNA data..."
	cancer_type = directory.split('/')[-1]

	matrix_file = os.path.join(directory, cancer_type +'.data.txt.matrix.npy')
	print 'Reading %s' %(matrix_file)
	data = np.load(open(matrix_file))
	patients = []
	genes = []

	gene_file = os.path.join(directory, cancer_type +'.data.txt.genes.csv')
	print 'Reading %s' %(gene_file)
	with open(gene_file, 'rb') as csvfile:
		reader = csv.reader(csvfile)
		count = 0
		for row in reader:
			for item in row:
				genes.append(item)
				count += 1

	patient_file = os.path.join(directory, cancer_type + '.data.txt.patients.csv')
	print 'Reading %s' %(patient_file)
	with open(patient_file, 'rb') as f:
		reader = csv.reader(f)
		patient_count = 0
		for row in reader:
			for item in row:
				patient_count += 1
				patients.append('-'.join(item.strip().split('-')[:3]))

	# Patients are columns
	assert patient_count == len(data[0])
	# Genes are rows
	assert count == data.shape[0]

	clinical_file = os.path.join(directory, cancer_type +  '.clinical/' + cancer_type + '.clin.merged.picked.txt.saved.p')
	print 'Reading %s' %(clinical_file)
	clinical = pickle.load(open(clinical_file))

	print 'All saved RNA files read'

	races = []
	patient_list = list(patients)
	indices_to_delete = set()
	c = 0
	for j, patient in enumerate(patient_list):
		if patient in clinical:
			info = clinical[patient]
			if info['race'] == 'NA':
				indices_to_delete.add(j)
				c +=1
			else:
				races.append(info['race'])
				c+=1
		else:
			indices_to_delete.add(j)
			c +=1

	assert(c == len(patient_list))

	# I delete patients that do not have clinical information. 
	# I should mask the corresponding rows in the matrix.
	patient_list = [item for a, item in enumerate(patient_list) if a not in indices_to_delete]
	# patient_list is now the set of patients to include in the study.

	if indices_to_delete:
		# Create mask to exclude patient data at indices to delete in our data matrix.
		mask = np.zeros_like(data)
		idx = np.array(list(indices_to_delete))
		mask[:,idx] = 1
		mask = np.ma.make_mask(mask)
		assert(np.ma.is_mask(mask))

		# Mask data matrix
		data = np.ma.masked_array(data, mask)
		assert(np.ma.is_masked(data))

	# Get the races of the patientisfo by reading clinical data.
	classes = races

	# Take the transpose to get into the format (n_samples, n_features)
	matrix = np.ma.transpose(data)
	
	# After transpose, num columns should match num genes and num rows shuld match number of patients.
	assert(matrix.shape[1] == len(genes))
	# Get the names of the genes (feature labels)
	labels = genes	
	return matrix, classes, labels
