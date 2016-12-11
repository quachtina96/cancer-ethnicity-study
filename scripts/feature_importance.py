'''fileoverview This file handles the analysis of feature importances that 
result from training classifiers on data from TCGA (The Cancer Genome Atlas) 
on GDC (Genomic Data Commons).

@author Tina Quach (quacht@mit.edu)
'''

import numpy as np 
from collections import defaultdict
import os
import subprocess
import sys
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt

class FeatureImportances:
	def __init__(self):
		self.numpy_array = None
		self.numpy_file = None
		self.feature_labels = None

	def set_importances(self, importances):
		self.numpy_array = importances
		self.numpy_file = None
	
	def set_feature_labels(self, labels):
		self.feature_labels = labels
	
	def save_npy(self, feature_importances_file, opt_feature_labels_file=None):
		'''Save feature importances (and labels) to numpy files'''
		np.save(feature_importances_file, self.numpy_array)
		if opt_feature_labels_file:
			np.save(opt_feature_labels_file, self.feature_labels)


	def load_npy(self, feature_importances_file, opt_feature_labels_file=None):
		# Expects the file inputs to be .npy files, loadable by numpy.
		self.numpy_array = np.load(feature_importances_file)
		self.numpy_file = feature_importances_file
		if opt_feature_labels_file:
			self.feature_labels = np.load(opt_feature_labels_file)
	
	def get_mean(self):
		return np.mean(self.numpy_array)

	def get_stdev(self):
		return np.std(self.numpy_array)

	def get_stats(self):
		stats_dict = defaultdict(float)
		stats_dict['Standard Deviation'] = self.get_stdev()
		stats_dict['Max'] = max(self.numpy_array)
		stats_dict['Min'] = min(self.numpy_array)
		stats_dict['Median'] = np.median(self.numpy_array)
		stats_dict['Mean'] = self.get_mean()
		return stats_dict
	
	def get_most_important_features(self, n_features, sorted_greatest_to_least=True):
		if not self.feature_labels:
			print "Cannot return features. First set_features"
			return
		else:
			return [self.feature_labels[indices] for indices in get_most_important_features_indices(n_features, sorted_greatest_to_least)]			

	def get_most_important_features_indices(self, n_features, sorted_greatest_to_least=True):
		# sorted least_to_greatest
		indices = self.numpy_array.argsort()[-n_features:]
		if sorted_greatest_to_least:
			# reverse the sort
			indices = indices[::-1]
		return indices
	
	def get_least_important_features_indices(self, n_features, sorted_greatest_to_least=True):
		# sorted least_to_greatest
		indices = self.numpy_array.argsort()[:n_features]
		if sorted_greatest_to_least:
			# reverse the sort
			indices = indices[::-1]
		return indices
			
	def pretty_print_map(self, n_features):
		sorted_important_indices = self.get_most_important_features_indices(n_features)
		for index in sorted_important_indices:
			print ('%f %s') %(self.numpy_array[index], self.feature_labels[index])

	def get_feature_importance_map(self, n_features):
		feature_importance_map = dict()
		indices = self.get_most_important_features_indices(n_features)
		for index in indices:
			feature = self.feature_labels[index]
			importance = self.numpy_array[index]
			feature_importance_map[feature] = importance
		return feature_importance_map

	def save_plot(self, n_features, plotfile):
		# Create barplot of the SNPS and the feature importances
		y_pos = np.arange(len(self.feature_labels))
		importances = self.numpy_array[self.get_most_important_features_indices()]

		plt.bar(y_pos, importances, align='center', alpha=0.5)
		plt.xticks(y_pos, self.feature_labels)
		plt.ylabel('Feature Importances')
		plt.title('Feature Importances of top ' + n_features +' SNPs')
		plt.savefig(plotfile)

if __name__ == '__main__':
	# parse command-line arguments
	if len(sys.argv) < 1:
		print "you must call program as:  "
		print "   python feature_importance.py <classifier_file> <labels npy file>"
		sys.exit(1)
	classifier_file = sys.argv[1]	
	labels_file = sys.arv[2]

	#Load classifier and labels
	classifier = pickle.load(open(classifer_file, "rb" ))
	labels = np.load(labels_file)	
	
	fi = FeatureImportances()
	fi.set_importances(classifier.feature_importances_)
	fi.set__feature_labels(labels)
	fi_stats = fi.get_stats()

	# Print out analysis of tree depths
	print'Feature Importances'
	for stat in fi_stats:
		print "%s: %f" %(stat, stats[stat])
	
	# Analyze feature importance to get the best
	print "The " + n_features + "Most Important SNPs"
	n_features = 10
	feature_to_importance = fi.get_feature_importance_map(10)
	fi.pretty_print_map(10)




