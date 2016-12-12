import csv
import pandas as pd
import numpy as np
import pickle
import statsmodels.api as sm
from statsmodels.formula.api import ols
import sys
from statsmodels.sandbox.stats.multicomp import fdrcorrection0, multipletests

def load_data(filename, cancer):
	data = np.load(open(filename +'.data.txt.matrix.npy'))
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

	assert patient_count == len(data[0])
	assert count == data.shape[0]

	clinical = pickle.load(open(filename + '.clinical/' + cancer + '.clin.merged.picked.txt.saved.p'))

	print 'All saved files read'

	index = list(patients)
	patient_data = {
			'age': [],
			'gender': [],
			'race': [],
			#'stage': []
	}

	indices_to_delete = set()
	for j, patient in enumerate(index):
		if patient in clinical:
			info = clinical[patient]
			if 'NA' in info.values():
				indices_to_delete.add(int(j))
			else:
				patient_data['age'].append(int(info['years_to_birth']))
				patient_data['gender'].append(info['gender'])
				patient_data['race'].append(info['race'])
				#patient_data['stage'].append(info['pathologic_stage'])
		else:
			indices_to_delete.add(int(j))

	index = [item for a, item in enumerate(index) if a not in indices_to_delete]
	print 'NUmber of cases deleted: ' + str(len(indices_to_delete))
	print 'Number of cases remaining: ' + str(len(index))

	assert len(index) == len(patient_data['gender'])

	master_df = pd.DataFrame(patient_data, index=index)
	print 'Master dataframe created'

	return data, patients, genes, master_df, clinical, indices_to_delete

def run_all(filename, cancer):
	data, patients, genes, master_df, clinical, indices_to_delete = load_data(filename, cancer)

	outfile = open(filename + '_pvals.tsv', 'w')
	for i, gene in enumerate(genes):
		expression_data = [item for a, item in enumerate(list(data[i])) if a not in indices_to_delete]
		df = master_df.copy()
		df['expression'] = pd.Series(expression_data, index=master_df.index)
		assert len(df['expression']) == len(df['gender']) == len(df['race']) == len(df['age'])
		# Run anova
		pval = anova(df)
		outfile.write(str(gene) + "\t" + str(pval) + "\n")

	print 'Completed ANOVA on all genes'
	print 'Output available at ' + filename + '_pvals.tsv'

def anova(data):
	#mod = ols('expression ~ age + gender + C(race) + C(stage)', data=data).fit()
	mod = ols('expression ~ age + gender + C(race)', data=data).fit()
	aov_table = sm.stats.anova_lm(mod)
	return float(aov_table.loc['C(race)']['PR(>F)'])

def main():

    """Checks if we have the right number of command line arguments
       and reads them in"""
    if len(sys.argv) < 1:
        print "you must call program as: python ./glm.py <datafile> <cancertype>"
        sys.exit(1)

    filename = sys.argv[1]
    cancer_type = sys.argv[2]
    run_all(filename, cancer_type)

if __name__ == "__main__":
    main()

#run_all('../data/BRCA/BRCA', 'BRCA')

