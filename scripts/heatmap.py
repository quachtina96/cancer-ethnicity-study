import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
import seaborn as sns
import csv
import sys

def plot_heatmap(filename, cancer):
	data = np.load(open(filename + '.data.txt.matrix.npy'))
	patients = []
	genes = []

	with open(filename +'_finalsnps.txt', 'rb') as genefile:
		for line in genefile:
			parts = line.strip().split('\t')
			tup = (int(parts[0]), parts[1], float(parts[2]))

			genes.append(tup)

	with open(filename + '.data.txt.patients.csv', 'rb') as f:
		reader = csv.reader(f)
		patient_count = 0
		for row in reader:
			for item in row:
				patient_count += 1
				patients.append('-'.join(item.strip().split('-')[:3]))
	genes = genes[:10]
	clinical = pickle.load(open(filename + '.clinical/' + cancer + '.clin.merged.picked.txt.saved.p'))
	reordered = []
	baa = []
	na = []
	ai = []
	white = []
	asian = []

	for i, patient in enumerate(patients):
		rae = 'NA'
		if patient in clinical:
			race = clinical[patient]['race'].strip()
		if race == 'white':
			white.append((i, patient))
		elif race == 'black or african american':
			baa.append((i, patient))
		elif race == 'american indian or alaska native':
			ai.append((i, patient))
		elif race == 'asian':
			asian.append((i, patient))
		else:
			na.append((i, patient))
	reordered = white + asian + baa + ai + na

	a = []
	for i, patient in enumerate(patients):
		for gene in genes:
			a.append(i)

	b = []
	for patient in patients:
		for ind, gene, pval in genes:
			b.append(gene)

	c = []
	for i, patient in reordered:
		for ind, gene, pval in genes:
			c.append(data[ind][i])

	matrix = {
		'patient': a,
		'gene': b,
		'expression': c
	}


	df = pd.DataFrame(matrix, columns=['patient', 'gene', 'expression'])
	df = pd.pivot_table(df,values='expression',index='gene',columns='patient')

	ax = sns.heatmap(df, xticklabels = 200, cmap='Reds')

	ax.set_title(cancer + " Gene Expression Heatmap")
	ax.vlines([len(white), len(white) + len(asian),
		len(white) + len(asian) + len(baa),
		len(white)+len(asian)+len(baa)+len(ai), len(white)+len(asian)+len(baa)+len(ai)+len(na)], *ax.get_ylim())
	sns.plt.yticks(rotation=0) 
	sns.plt.xticks(rotation=270)
	sns.plt.savefig(filename + "_heatmap.png")

def main():

    """Checks if we have the right number of command line arguments
       and reads them in"""
    if len(sys.argv) < 1:
        print "you must call program as: python ./glm.py <datafile> <cancertype>"
        sys.exit(1)

    filename = sys.argv[1]
    cancer_type = sys.argv[2]
    plot_heatmap(filename, cancer_type)

if __name__ == "__main__":
    main()
