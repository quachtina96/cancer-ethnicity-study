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
import pickle
from feature_importance import FeatureImportances

class RandomForest:
	def __init__(self):
		self.classifier = None
		self.tree_depths = None

	def fit(self, data, classes):
		self.classifier = RandomForestClassifier(n_estimators=1000, oob_score=True, n_jobs=-1, verbose=1)
		self.classifier = self.classifier.fit(data, classes)
		self.tree_depths = self.get_tree_depths()
		return self.classifier

	def save(self, file):
		if self.classifier:
			p_file = open('/'.join(matrix_path.split('/')[:-1]) +'classifier.RF.p','w')
			p = pickle.dump(self.classifier, p_file, protocol=2)
		else:
			print 'Could not save classifier. Empty file.'

	def load(self, pickle_file):
		try:
			self.classifier = pickle.load(open(pickle_file, 'r'))
			self.tree_depths = self.get_tree_depths()
		except:
			print 'Could not load classifier'

	def get_tree_depths(self):
		return [estimator.tree_.max_depth for estimator in self.classifier.estimators_]

	def iterfit(self, data, classes):
		pass

	def quick_stats(self):
		print "Classes: %s" %(str(self.classifier.classes_))
		print "OOB Score: %f" %(self.classifier.oob_score_)
		print "Number of Outputs: %d" %(self.classifier.n_outputs_)

	def analyze_tree_depths(self):
		tree_depths = self.get_tree_depths()
		tree_depth_stats = defaultdict(float)
		tree_depth_stats['Max'] = max(tree_depths)
		tree_depth_stats['Min'] = min(tree_depths)
		tree_depth_stats['Median'] = min(tree_depths)
		tree_depth_stats['Mean'] = min(tree_depths)
		return tree_depth_stats

		
if __name__ == '__main__':
	# parse command-line arguments
	if len(sys.argv) < 1:
		print "you must call program as:  "
		print "   python randomforest.py <matrix_dir>  <matrix_type> <n_features>"
		print "   or "
		print "   python randomforest.py -h"
		sys.exit(1)

	if (sys.argv[1] == '-h'):
		# Display help message
		print "python randomforest.py <matrix_dir>  <matrix_type> <n_features>"
		print "   matrix_dir: path to the directory containing the matrix to train on. "
		print "   matrix_type: type of matrix ('snp')"
	directory = sys.argv[1]
	matrix_type = sys.argv[2]
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
