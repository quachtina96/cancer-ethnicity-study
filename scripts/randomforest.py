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
import matplotlib.pyplot as plt

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
			if file.endswith("matrix.tsv"):
				print file
				matrix_list.append(os.path.join(root, file))
				print os.path.join(root, file)
	
	matrix_path = matrix_list[0]
	matrix_filename = matrix_path.split('/')[-1]
	
	snp_matrix = np.genfromtxt(matrix_path, delimiter='\t', dtype=None,names=True)
	labels = snp_matrix.dtype.names
	print labels
	print "Number of patients before filtering those without reported races: %d" %(snp_matrix.shape[0])
	
	# Filter out the snp_matrix such that only elements that have specified
	# races are included in the classification. At the same time convert 
	# snp_matrix from being a np.array of tuples to np.array of lists (which has
	# the correct shape for analysis).
	# if not line['Race'] try line[-2]
	filtered = np.array([list(line) for line in snp_matrix if line['Race'] != 'not reported'])
	print "Number of patients after filtering: %d" %(filtered.shape[0])

	snp_data = filtered[:,1:-2]
	print "Number of SNPs considered: %d" %(snp_data.shape[1])
	races = filtered[:,-2]	

	print 'Training RandomForestClassification on given data'
	#TODO: make this nonverbose
	random_forest = RandomForestClassifier(n_estimators=1000, oob_score=True, n_jobs=-1, verbose=1)
	fitted = random_forest.fit(snp_data,races)
	tree_depths = [estimator.tree_.max_depth for estimator in fitted.estimators_]

	fittedInfoDict = defaultdict(dict)

	print "RandomForestClassifier statistics:"
	print "Classes: %s" %(str(fitted.classes_))
	print "OOB Score: %f" %(fitted.oob_score_)
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
	sorted_important_indices = 
	for indices in most_important_indices[::-1]:
		print ('%f %s') %(fitted.feature_importances_[indices], labels[indices+1])
	
	# Create barplot of the SNPS and the feature importances
	objects = labels[[index + 1 for index in most_important_indices[::-1]]]
	y_pos = fitted.feature_importances_[most_important_indiceis[]

performance = [10,8,6,4,2,1]
 
plt.bar(y_pos, performance, align='center', alpha=0.5)
plt.xticks(y_pos, objects)
plt.ylabel('Usage')
plt.title('Programming language usage')
 
plt.show()]
performance = [10,8,6,4,2,1]
 
plt.bar(y_pos, performance, align='center', alpha=0.5)
plt.xticks(y_pos, objects)
plt.ylabel('Usage')
plt.title('Programming language usage')
 
plt.show()
