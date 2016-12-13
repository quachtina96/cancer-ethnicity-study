import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import pickle
import csv

filename = '../data/BRCA/BRCA'
cancer = 'BRCA'
data = np.load(open(filename +'.data.txt.matrix.npy')).transpose()
pca = PCA(n_components=2)
reduced = pca.fit_transform(data)
clinical = pickle.load(open(filename + '.clinical/' + cancer + '.clin.merged.picked.txt.saved.p'))
patients = []
with open(filename + '.data.txt.patients.csv', 'rb') as f:
	reader = csv.reader(f)
	patient_count = 0
	for row in reader:
		for item in row:
			patient_count += 1
			patients.append('-'.join(item.strip().split('-')[:3]))

plt.figure()
for i, item in enumerate(reduced):
	color = 'b'
	if patients[i] in clinical:
		if clinical[patients[i]]['race'].strip() == 'white':
			color = 'r'
		elif clinical[patients[i]]['race'].strip() == 'asian':
			color = 'm'
		elif clinical[patients[i]]['race'].strip() == 'black or african american':
			color = 'y'
	plt.plot(item[0], item[1], 'ro', color=color)
# plt.scatter(reduced[:,0], reduced[:,1])
plt.xlabel('1st eigenvector')
plt.ylabel('2nd eigenvector')
plt.show()


