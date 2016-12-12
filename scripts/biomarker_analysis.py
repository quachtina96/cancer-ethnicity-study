'''
@fileoverview This file analyzes the biological significance of biomarkers
identified through feature selection on data from TCGA (The Cancer Genome Atlas)

@author Tina Quach (quacht@mit.edu)
'''
import mini
from feature_importance import FeatureImportances
from imblearn.ensemble import EasyEnsemble

def readSNP(matrix_path):
	# Read SNP Matrix
	snp_matrix = np.genfromtxt(matrix_path, delimiter='\t', dtype=None,names=True)
	labels = snp_matrix.dtype.names[1:-2]
	print "Number of patients before filtering those without reported races: %d" %(snp_matrix.shape[0])

	# Filter out the snp_matrix such that only elements that have specified
	# races are included in the classification. At the same time convert 
	# snp_matrix from being a np.array of tuples to np.array of lists (which has
	# the correct shape for analysis).
	# if not line['Race'] try line[-2]
	filtered = np.array([list(line) for line in snp_matrix if line['Race'] != 'not reported'])
	print "Number of patients after filtering: %d" %(filtered.shape[0])

	data = filtered[:,1:-2]

	print "Number of SNPs considered: %d" %(data.shape[1])
	races = filtered[:,-2]
	return (data, races)

def readSNPnpy(matrix_path):
	# Read SNP Matrix
	snp_matrix = np.load(matrix_path)
	# Filter out the snp_matrix such that only elements that have specified
	# races are included in the classification. At the same time convert 
	# snp_matrix from being a np.array of tuples to np.array of lists (which has
	# the correct shape for analysis).
	# if not line['Race'] try line[-2]
	filtered = np.array([list(line) for line in snp_matrix if line['Race'] != 'not reported'])
	print "Number of patients after filtering: %d" %(filtered.shape[0])

	data = filtered[:,1:-2]

	print "Number of SNPs considered: %d" %(data.shape[1])
	races = filtered[:,-2]
	return (data, races)

'''
Find all files in the directory that have feature importance in its name.
Get that list of files.
For every file in this list,
Parse out CANCER TYPE from the file name?
load and cast to a feature importances object.
# '''

# cancer_types = mini.mini_labels
# importance_info_list = mini.mini_dict_list()

# for cancer_type in importance_info_list:
# 	# load appropriate matrix
# 	cancer_type[]

# # We are interested in across all cancer types.
# # For every selected feature

# # What proportion of each race or ethnicity has the specific trait

# # What is the distrubution of races for that particular S

# Generate the 3 dimensional resampled data where first dimension refers to 
# the number of subsets generated by the EasyEnsemble.

if __name__ == '__main__':
	# parse command-line arguments
	if len(sys.argv) < 1:
		print "you must call program as:  "
		print "   python randomforest.py <matrix_path>"
		sys.exit(1)

	matrix_path = sys.argv[1]
	if matrix_path.split('.')[-1] == 'npy':
		data, classes = readSNPnpy(matrix_path)
	else:
		data, classes = readSNP(matrix_path)
	
	
		
	n_subsets = 10
	ee = EasyEnsemble(n_subsets)
	ee.fit(data, classes)
	print ee.stats_c_
