import numpy as np
import csv

# Each column is a patient ID
# First column is a gene
# Up to line 31 are unknown genes

filename = '../data/BRCA/BRCA.txt'

patients = None
index = 0
end_range = -1
genes = []
num_skip = None

with open(filename) as f:
	for line in f:
		if index == 0:
			patients = line.split("\t")[1:]
			end_range = len(patients) + 1

		if index > 1:
			parts = line.split("\t")
			gene_id = parts[0]
			if gene_id[0] != "?":
				if num_skip == None:
					num_skip = index
				genes.append(gene_id)
		index += 1

data = np.loadtxt(filename, delimiter='\t', skiprows=num_skip, usecols=range(1,end_range))

out = open('BRCA_RNA_patients.csv', 'wb')
wr = csv.writer(out, quoting=csv.QUOTE_ALL)
wr.writerow(patients)

# Transpose matrix (row is patient, column is gene)
data = data.transpose()

assert len(genes) == data.shape[1]
assert len(patients) == data.shape[0]

np.save('BRCA_RNA_matrix.out', data)