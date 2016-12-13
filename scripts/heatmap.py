import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import csv

# You want to choose top 20 genes
# Get the expression data for those 20 genes
# Plot that shiz


filename = '../data/BRCA/BRCA'
data = np.load(open('../data/BRCA/BRCA.data.txt.matrix.npy'))
patients = []
genes = []

# month = gene
# year = patient
# passengers = level

with open(filename +'.data.txt.genes.csv', 'rb') as csvfile:
	reader = csv.reader(csvfile)
	count = 0
	for row in reader:
		for item in row:
			genes.append(item)
			count += 1

with open(filename + '.data.txt.patients.csv', 'rb') as f:
	reader = csv.reader(f)
	patient_count = 0
	for row in reader:
		for item in row:
			patient_count += 1
			patients.append('-'.join(item.strip().split('-')[:3]))

patients = patients[:60]
genes = genes[:20]

a = []
for patient in patients:
	for gene in genes:
		a.append(patient)
b = []
for patient in patients:
	for gene in genes:
		b.append(gene)

c = []
for i, patient in enumerate(patients):
	for j, gene in enumerate(genes):
		c.append(data[j][i])

matrix = {
	'patient': a,
	'gene': b,
	'expression': c
}


df = pd.DataFrame(matrix, columns=['patient', 'gene', 'expression'])

df = df.pivot("gene", "patient", "expression")
ax = sns.heatmap(df)
sns.plt.yticks(rotation=0) 
sns.plt.xticks(rotation=270)
sns.plt.show()