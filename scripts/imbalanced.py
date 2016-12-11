'''fileoverview This file handles iterative random forest analysis on cancer
data from TCGA (The Cancer Genome Atlas) on GDC (Genomic Data Commons).

@author Tina Quach (quacht@mit.edu)
'''

# import numpy as np 
# from randomforest import RandomForest
# from collections import defaultdict
# from imblearn.ensemble import EasyEnsemble
# from feature_importance import FeatureImportances

import matplotlib.pyplot as plt
import numpy as np
from sklearn.datasets import make_classification
from sklearn.decomposition import PCA
from imblearn.ensemble import EasyEnsemble

# Apply Easy Ensemble
ee = EasyEnsemble()
X_resampled, y_resampled = ee.fit_sample(X, y)

if __name__ == '__main__':
	directory = sys.argv[1]
	matrix_type = sys.argv[2]
	if sys.argv[3]:
		n_features = sys.argv[3]

	print "Conducting Random Forest analysis..."
	for root, dirs, files in os.walk(directory):
		for file in files:
			if file.endswith("matrix.tsv"):
				print file
				matrix_path = os.path.join(root, file)
	
	matrix_filename = matrix_path.split('/')[-1]
	
	if matrix_type == 'snp':
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
		classes = races
		
		matrix = snp_matrix	


	print 'Training RandomForestClassifier on given data...'
	rf = RandomForest()
	rf.fit(data, classes)

	print 'Saving the classifier...'
	p_file = open('/'.join(matrix_path.split('/')[:-1]) +'classifier.RF.p','w')
	rf.save(p_file)

	# Analyze Classifier
	rf.quick_stats()
	tree_depth_stats = rf.analyze_tree_depths()
	
	# Write all feature importances and original data matrix to a file
	np.save(matrix_path + '_RF_feature_importance', rf.classifier.feature_importances_)
	np.save(matrix_path + '_snps', labels)
	np.save(matrix_path, matrix)

	if (n_features)
		fi = FeatureImportances()
		fi.set_importances(rf.classifier.feature_importances_)
		fi.set_feature_labels(labels)
		fi_stats = fi.get_stats()

		# Print out analysis of tree depths
		infoDict = {'Tree_Depth': tree_depth_stats, 'Feature Importances': fi_stats}
		for attribute in infoDict:
			print attribute
			stats = infoDict[attribute]
			for stat in stats:
				print "%s: %f" %(stat, stats[stat])
		
		# Analyze feature importance to get the best
		print "The " + str(n_features) + " Most Important SNPs"
		feature_to_importance = fi.get_feature_importance_map(n_features)
		fi.pretty_print_map(n_features)