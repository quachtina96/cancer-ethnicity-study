'''fileoverview This file handles the analysis of feature importances that result from training classifiers on data from TCGA (The Cancer Genome Atlas) on GDC (Genomic Data Commons).

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
	def __init__(self, feature_importances_file, opt_feature_labels_file=None):
		# Expects the file inputs to be .npy files, loadable by numpy.
		self.numpy_array = np.load(feature_importances_file)
		self.numpy_file = feature_importances_file
		if opt_feature_labels_file:
			self.feature_labels = np.load(opt_feature_labels_file)
	
	def get_mean(self):
		return np.mean(self.numpy_array)

	def get_stdev(self):
		return np.std(self.numpy_array)
	
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
			
	
	def get_feature_importance_map(self, n_features):
		feature_importance_map = dict()
		indices = get_most_important_features_indices(n_features)
		for index in indices:
			feature = self.feature_labels[index]
			importance = self.numpy_array[index]
			feature_importance_map[feature] = importance
		return feature_importance_map

