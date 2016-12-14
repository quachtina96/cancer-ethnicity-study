import numpy as np 
import csv

for cancer in ['BRCA', 'LGG', 'GBM', 'LUAD', 'LUSC']:
	print cancer
	filename = '../data/rna_data/' + cancer + '/'+ cancer
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
	print count
	print ""