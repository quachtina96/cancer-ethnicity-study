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

class RandomForest:
	def __init__(self):
		self.classifier = None
		self.tree_depths = None

	def fit(self, data, classes):
		if not self.classifier:
			self.classifier = RandomForestClassifier(n_estimators=1000, oob_score=True, n_jobs=-1, verbose=1)
		self.classifier = self.classifier.fit(data, classes)
		self.get_tree_depths = self.get_tree_depths()
		return self.classifier

	def save(self, file)
		if self.classifier:
			p_file = open('/'.join(matrix_path.split('/')[:-1]) +'classifier.RF.p','w')
			p = pickle.dumps(self.classifier, p_file)
		else:
			print 'Could not save classifier. Empty file.'

	def load(self, pickle_file):
		try:
			p = pickle.dumps(fitted, p_file)
			self.classifier = pickle.loads(open(pickle_file, 'r'))
			self.get_tree_depths = self.get_tree_depths()
		except:
			print 'Could not load classifier'

	def get_tree_depths(self):
		return [estimator.tree_.max_depth for estimator in self.classifier.estimators_]

	def iterfit(self, data, classes):
		pass

	def analyze(self):
		print "Classes: %s" %(str(self.classifier.classes_))
		print "OOB Score: %f" %(self.classifier.oob_score_)
		print "Number of Outputs: %d" %(self.classifier.n_outputs_)

		fitted_attributes = {'Feature Importance': self.classifier.feature_importances_,'Tree Depths':self.get_tree_depths()}
		fittedInfoDict = defaultdict(dict)
		for interest in fitted_attributes:
			fittedInfoDict[interest]['Max'] = max(fitted_attributes[interest])
			fittedInfoDict[interest]['Min'] = min(fitted_attributes[interest])
			fittedInfoDict[interest]['Median'] = min(fitted_attributes[interest])
			fittedInfoDict[interest]['Mean'] = min(fitted_attributes[interest])

		for attribute in fittedInfoDict:
			attributeInfo = fittedInfoDict[attribute]
			for stat in attributeInfo:
				print "%s: %f" %(stat, attributeInfo[stat])
		return fittedInfoDict

if __name__ == '__main__':
	# parse command-line arguments
	if len(sys.argv) < 1:
		print "you must call program as:  "
		print "   python randomforest.py <matrix_dir>  <matrix_type>"
		sys.exit(1)
	directory = sys.argv[1]
	matrix_type = sys.argv[2]

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
		labels = snp_matrix.dtype.names
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

	print 'Training RandomForestClassifier on given data...'
	rf = RandomForest(data,races)

	print 'Saving the classifier...'
	p_file = open('/'.join(matrix_path.split('/')[:-1]) +'classifier.RF.p','w')
	rf.save(p_file)

	# Analyze Classifier
	fittedInfoDict = rf.analyze()
	for attribute in fittedInfoDict:
		print attribute
		attributeInfo = fittedInfoDict[attribute]
		for stat in attributeInfo:
			print "%s: %f" %(stat, attributeInfo[stat])
	
	# Analyze feature importance to get the best
	print "The Most Important SNPs"
	snp_labels = []
	most_important_indices = fitted.feature_importances_.argsort()[-10:]
	sorted_important_indices = most_important_indices[::-1]
	for indices in sorted_important_indices:
		print ('%f %s') %(fitted.feature_importances_[indices], labels[indices+1])
		snp_labels.append(labels[indices+1])

	# Write all feature importances and SNPs to a file
	np.save(matrix_path + '_RF_feature_importance', fitted.feature_importances_)
	np.save(matrix_path + '_snps', labels[1:-2])
	np.save(matrix_path, snp_matrix)

	# Create barplot of the SNPS and the feature importances
	# y_pos = np.arange(len(snp_labels))
	# importances = fitted.feature_importances_[sorted_important_indices]

	# plt.bar(y_pos, importances, align='center', alpha=0.5)
	# plt.xticks(y_pos, snp_labels)
	# plt.ylabel('Feature Importances')
	# plt.title('Feature Importances of top 10 SNPs')
	# plt.savefig(matrix_path +'_barplot.png')

