'''
@fileoverview This file handles the random forest analysis on cancer
data from TCGA (The Cancer Genome Atlas) on GDC (Genomic Data Commons).

@author Tina Quach (quacht@mit.edu)
'''

import numpy as np 
from collections import defaultdict
import os
import sys
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
import pickle
from feature_importance import FeatureImportances
from imblearn.ensemble import EasyEnsemble
from randomforest import RandomForest

N_ESTIMATORS = 100
def run_random_forest_analysis(data, classes, matrix_dir, analysisID):	
	print 'Training RandomForestClassifier on given data...'
	rf = RandomForest(N_ESTIMATORS)
	rf.fit(data, classes)

	print 'Saving the classifier...'
	classifier_file = str(analysisID) + 'classifier.RF.p'

	p_file = os.path.join(matrix_dir,classifier_file)
	rf.save(p_file)

	# Analyze Classifier
	if rf.classifier.oob_score_ == 0:
		print rf.classifier.get_params()
	rf.quick_stats()
	tree_depth_stats = rf.analyze_tree_depths()

	# Write all feature importances and original data matrix to a file
	np.save(matrix_path + '_RF_feature_importance', rf.classifier.feature_importances_)
	print "Classifier saved to %s" %(p_file)
	print "Feature Importances saved to %s" %(matrix_path + '_RF_feature_importance.npy')
	return rf

if __name__ == '__main__':
	# parse command-line arguments
	if len(sys.argv) < 2:
		print "you must call program as:  "
		print "   python randomforest.py <matrix_dir>  <matrix_type> <opt_n_features>"
		print "where"
		print "   matrix_dir: path to the directory containing the matrix to train on. "
		print "   matrix_type: type of matrix ('snp')"
		print "   opt_n_features: number of features to identify"
		sys.exit(1)

	directory = sys.argv[1]
	matrix_type = sys.argv[2]
	if sys.argv[3]:
		n_features = int(sys.argv[3])

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


	# TODO: handle multiple matrix types here.
	np.save(matrix_path + '_snps', labels)
	np.save(matrix_path, matrix)

	print 'Sampling imbalanced data using EasyEnsemble to generate balanced \
	samples...'

	# Apply Easy Ensemble
	n_subsets = 10
	ee = EasyEnsemble(n_subsets=n_subsets)

	# Generate the 3 dimensional resampled data where first dimension refers to 
	# the number of subsets generated by the EasyEnsemble.
	X_resampled, y_resampled = ee.fit_sample(data, classes)
	
	for i in xrange(n_subsets):
		analysisID = str(i)
		data = X_resampled[i]
		classes = y_resampled[i]
		rf = run_random_forest_analysis(data, classes, directory, analysisID)

		if n_features:
			fi = FeatureImportances()
			fi.set_importances(rf.classifier.feature_importances_)
			fi.set_feature_labels(labels)
			fi_stats = fi.get_stats()

			# Print out analysis of tree depths
			tree_depth_stats = rf.analyze_tree_depths()
			infoDict = {'Tree_Depth': tree_depth_stats, 'Feature Importances': fi_stats}
			for attribute in infoDict:
				print attribute
				stats = infoDict[attribute]
				for stat in stats:
					print "%s: %f" %(stat, stats[stat])
			
			# Analyze feature importance to get the best
			print "The %d Most Important SNPs" %(n_features)
			feature_to_importance = fi.get_feature_importance_map(n_features)
			fi.pretty_print_map(n_features)
