import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import csv
import pickle
import os 

class SNP:
	def __init__(self, label):
		self.label = label 
		temp = label.split('_')
		self.gene = label[0]
		self.chr = label[1]
		self.start_pos = label[2]

def read_mini(mini):
	scores = []
	features = []
	featuress = mini.strip().split('\n')
	for feature in featuress:
		splitted = feature.split(' ')
		score = splitted[0]
		name = splitted[1]
		scores.append(float(score))
		features.append(name)
	return (scores,features)

# Get selected SNPs:
minidir = '../results/randomforest'
minidir = '../results/rna_plot'
cancer_types = ['BRCA', 'GBM', 'LGG', 'LUAD', 'LUSC']

for cancer in cancer_types:
	minifile = cancer + '.rf.selected.txt'
	with open(os.path.join(minidir,minifile), 'r') as content_file:
	    content = content_file.read()
	scores, features = read_mini(content)
	selected = features[:20]
	# You want to choose top 20 genes
	# Plot that shiz
	cancer = 'BRCA'
	filename = '../data/rna_data/' + cancer
	data = np.load(open(filename + '.data.txt.matrix.npy'))
	patients = []
	genes = []


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
		for gene in selected:
			a.append(i)

	b = []
	for patient in patients:
		for gene in selected:
			b.append(gene)

	# c = []
	# for i, patient in enumerate(patients):
	# 	for ind, gene, pval in genes:
	# 		c.append(data[ind][i])

	c = []
	for i, patient in reordered:
		for gene in selected:
			ind = list(genes).index(gene)
			c.append(data[ind][i])

	matrix = {
		'patient': a,
		'gene': b,
		'expression': c
	}


	df = pd.DataFrame(matrix, columns=['patient', 'gene', 'expression'])
	df = pd.pivot_table(df,values='expression',index='gene',columns='patient')

	ax = sns.heatmap(df, xticklabels = 200)

	sns.plt.yticks(rotation=0) 
	sns.plt.xticks(rotation=270)
	sns.plt.tight_layout()
	sns.plt.savefig(os.path.join(minidir,minifile+'heatmap.png'))
