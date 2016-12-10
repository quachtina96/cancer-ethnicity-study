'''
@fileoverview This file handles the random forest analysis on cancer
data from TCGA (The Cancer Genome Atlas) on GDC (Genomic Data Commons).

@author Tina Quach (quacht@mit.edu)
'''

import numpy as np 
from collections import defaultdict
import os
import subprocess
import sys
from sklearn.ensemble import RandomForestClassifier

if __name__ == '__main__':
	# parse command-line arguments
	if len(sys.argv) < 1:
		print "you must call program as:  "
		print "   python randomforest.py <matrix_dir>"
		sys.exit(1)
	directory = sys.argv[1]

	matrix_list = []
	print "List of Matrix TSV files using for Random Forest analysis..."
	for root, dirs, files in os.walk(directory):
		for file in files:
			if file.endswith(".tsv"):
				print file
				matrix_list.append(os.path.join(root, file))
				print os.path.join(root, file)
	
	matrix_path = matrix_list[0]
	matrix_filename = matrix_path.split('/')[-1]
	
	snp_matrix = np.genfromtxt(matrix_path, delimiter='\t', dtype=None,names=True)
	print "shapes"
	print snp_matrix.shape
	# Filter out the snp_matrix such that only elements that have specified
	# races are included in the classification. At the same time convert 
	# snp_matrix from being a np.array of tuples to np.array of lists (which has
	# the correct shape for analysis).
	# if not line['Race'] try line[-2]
	filtered = np.array([list(line) if line['Race'] != 'not reported' for line in snp_matrix])
	print filtered.shape

	snp_data = filtered[:,1:-2]
	print snp_data.shape
	races = filtered[:,-2]
	print races.shape

	print 'Training RandomForestClassification on given data'
	#TODO: make this nonverbose
	random_forest = RandomForestClassifier(n_estimators=1000, oob_score=True, n_jobs=-1, verbose=1)
	fitted = random_forest.fit(snp_data,race)
	tree_depths = [estimator.tree_.max_depth for estimator in fitted.estimators_]

	fittedInfoDict = defaultdict(dict)

	print "RandomForestClassifier statistics:"
	print "Classes: %s" %(str(fitted.classes_))
	print "OOB Score: %d" %(fitted.oob_score_)
	print "Number of Outputs: %d" %(fitted.n_outputs_)

	classifier_interests = {'Feature Importance': fitted.feature_importances_,'Tree Depths':tree_depths}
	for interest in classifier_interests:
		classifier_interests[interest]
		fittedInfoDict['Max ' + interest] = max(classifier_interests[interest])
		fittedInfoDict['Min ' + interest] = min(classifier_interests[interest])
		fittedInfoDict['Median ' + interest] = min(classifier_interests[interest])
		fittedInfoDict['Mean ' + interest] = min(classifier_interests[interest])

	# Get Tree Depths 
	for attribute in fittedInfoDict:
		print "%s: %s" %(attribute, str(fittedInfoDict[attribute]))
	
	# Analyze feature importance to get the best
	print "The Most Important SNPS"
	most_important_indices = fitted.feature_importances_.argsort()[-10:]
	least_to_most = fitted.feature_importances_[most_important_indices]
	for snp in least_to_most[::-1]:
		print snp

