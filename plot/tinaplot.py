import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import csv

# You want to choose top 20 genes
# Get the expression data for those 20 genes
# Plot that shiz

cancer_type = 'BRCA'
#matrix_path = cancer_type + '_snp_matrix.tsv'
snp_matrix = np.load(matrix_path + '.npy')
snps = np.load('feature_labels.npy')
selected = ['AC0094876_chr2_161504316', 'INSIG2_chr2_118106806', 'LYST_chr1_235777293', 'AC0838993_chr2_87180029', 'ZDHHC17_chr12_76848315', 'RP13347D81_chrX_118938384', 'DOCK4_chr7_111989081', 'HLTF_chr3_149032250', 'LAPTM4B_chr8_97816170', 'CHRNA5_chr15_78589855', 'RP11307L32_chr9_6278929', 'ZNF510_chr9_96758941', 'NBPF13P_chr1_147044470', 'PTEN_chr10_87952218', 'ZNF223_chr19_44066431', 'AKT3_chr1_243552774', 'GRM1_chr6_146399052', 'STK39_chr2_168012666', 'RARS_chr5_168516842', 'MARC2_chr1_220762956', 'TSPYL1_chr6_116279518', 'TAS2R31_chr12_11030810', 'TCHH_chr1_152111278', 'CD300LG_chr17_43853893', 'EPHA1_chr7_143397954', 'MUC2_chr11_1098154', 'PCDH19_chrX_100350661', 'OXCT1_chr5_41731732']

filtered = np.array([list(line) for line in snp_matrix if line['Race'] != 'not reported'])

patients = filtered[:,0]
data = filtered[:,1:-2]
races = filtered[:,-2]
classes = races


# Construct the matrix to visualize
patients = patients[:60]
selected = selected[:20]

a = []
for patient in patients:
	for snp in selected:
		a.append(patient)

b = []
for patient in patients:
	for snp in selected:
		b.append(snp)

c = []
for i, patient in enumerate(patients):
	for snp in selected:
		j = snps.index(snp)
		c.append(data[i][j])

matrix = {
	'patient': a,
	'snp': b,
	'presence': c
}

# bottom, top left, top right
#columns = ['patient', 'snp', 'binary']
df = pd.DataFrame(matrix, columns=['patient', 'snp', 'presence'])

df = df.pivot("snp", "patient", "presence")
ax = sns.heatmap(df)
sns.plt.yticks(rotation=0) 
sns.plt.xticks(rotation=270)
sns.plt.show()
ax.savefig('seaborn_heatmap.png')


