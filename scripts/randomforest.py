'''
@fileoverview This file handles the random forest analysis on cancer
data from TCGA (The Cancer Genome Atlas) on GDC (Genomic Data Commons).

@author Tina Quach (quacht@mit.edu)
'''
import operator
import os
import subprocess
import sys
import numpy as np 
from collections import defaultdict
import pickle
from sklearn.ensemble import RandomForestClassifier
from feature_importance import FeatureImportances
import RNAreader
import pprint

class RandomForest:
	def __init__(self, n_estimators):
		self.classifier = None
		self.tree_depths = None
		self.n_estimators = n_estimators

	def fit(self, data, classes):
		self.classifier = RandomForestClassifier(n_estimators=self.n_estimators,n_jobs=-1, oob_score=True, verbose=1, class_weight="auto")
		self.classifier = self.classifier.fit(data, classes)
		self.tree_depths = self.get_tree_depths()
		return self.classifier

	def save(self, classifier_file):
		if self.classifier:
			try:
				p_file = open(classifier_file ,'w+')
				p = pickle.dump(self.classifier, p_file, protocol=2)
			except:
				'Could not save classifier.'
		else:
			print 'Could not save classifier.'

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

# Util functions
def renormalize(distr):
    normalization_constant = sum(distr.values())
    for key in distr.keys():
        distr[key] = float(distr[key])/normalization_constant
	
# Functions for comparing within a specific caner type.
def getClassProportions(data, classes_labels, index):
	'''What proportion of each race or ethnicity has the specific trait?'''

	# Get the distribution of a certain SNP's occurrence over the observed classes
	counts = getCountsOverClasses(data, classes_labels, index, False)
	proportions = defaultdict(float)
	for candidate_class in counts:
		class_member_count = list(classes_labels).count(candidate_class)
		proportions[candidate_class] = float(counts[candidate_class])/class_member_count
	return proportions

def getCountsOverClasses(data, class_labels, index, opt_normalized):
	'''What is the distribution of a certain SNP's occurrence over the observed
	races? This function gets the counts that you can normalize to get a 
	distribution'''
	feature_occurences = data[:, index]
	# print feature_occurences[:5]
	counts = defaultdict(float)	
	
	for i in xrange(data.shape[0]):
		if int(feature_occurences[i]):
			counts[class_labels[i]] += 1
	if opt_normalized:
		renormalize(counts)
	return counts
		
if __name__ == '__main__':
	# parse command-line arguments
	if len(sys.argv) < 2:
		print "you must call program as:  "
		print "   python randomforest.py <matrix_dir>  <matrix_type> <n_features>"
		print "an addition <n_iterations> to run random forest is optional"
	
		sys.exit(1)

	directory = sys.argv[1]
	matrix_type = sys.argv[2]
	n_features = int(sys.argv[3])
	n_iter = 1
	if len(sys.argv) == 4:
		n_iter = int(sys.argv[4])
	
	if matrix_type == 'snp':
		print "Conducting Random Forest analysis on SNP data..."
		for root, dirs, files in os.walk(directory):
			for file in files:
				if file.endswith("matrix.tsv"):
					print file
					matrix_path = os.path.join(root, file)
		
		matrix_filename = matrix_path.split('/')[-1]

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
		np.save(matrix_path, matrix)
	elif matrix_type == 'rna':

		matrix_path, data, classes, labels = RNAreader.read_rna(directory) 
		assert(np.ma.is_masked(data))
		data = np.array([line for line in data if not np.ma.is_masked(line)])
		classes = np.array(classes)
		matrix_path = '/home/cancer-ethnicity-study/results/randomforest/rna'
	else:
		print "Invalid matrix type. Must be 'snp' or 'rna'"
		sys.exit(1)

	print ""
	print 'Training RandomForestClassifier on data...'
	rf = RandomForest(1000)
	
	for i in xrange(n_iter):
		rf.fit(data, classes)

		print 'Saving the classifier...'
		try:
			p_file = '/'.join(matrix_path.split('/')) + 'classifier.RF.p'
			rf.save(p_file)
		except:
			print 'Could not save the classifier due to permissions'

		# Analyze Classifier
		rf.quick_stats()
		tree_depth_stats = rf.analyze_tree_depths()
		
		# Write all feature importances and original data matrix to a file
		feature_importance_file = matrix_path + '_RF_feature_importance_'+str(i)
		try:
			np.save(feature_importance_file, rf.classifier.feature_importances_)
		except:
			print 'Could not save feature importances'

	if matrix_type == 'snp':
		np.save(os.path.join(directory,'feature_labels'), labels)
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
		print ""
	
	# Analyze feature importance to get the best
	print "The %d Most Important Features" %(n_features)
	feature_to_importance = fi.get_feature_importance_map(n_features)
	fi.pretty_print_map(n_features)

	identified_features_info = defaultdict(defaultdict)
	for index in fi.get_most_important_features_indices(n_features):
		feature = labels[index] 
		proportions = getClassProportions(data, classes, index)
		counts = getCountsOverClasses(data, classes, index, False)
		identified_features_info[feature]['index'] = index
		identified_features_info[feature]['proportions'] = proportions
		identified_features_info[feature]['counts'] = counts

	# Filter the features based on proportions.
	feature_std_map = {}
	feature_list = identified_features_info.keys()
	for feature in feature_list:
		identified_features_info[feature]['std'] = np.std(identified_features_info[feature]['proportions'].values())
		if identified_features_info[feature]['std'] < .15:
			del identified_features_info[feature]
		else:
			feature_std_map[feature] = identified_features_info[feature]['std']

	# Sort the filtered features 
	sorted_by_std = sorted(feature_std_map.items(), key=operator.itemgetter(1))
	print ""
	print 'Identified %d features' %(len(feature_std_map.keys()))
	for item in reversed(sorted_by_std):
		feature = item[0]
		std = item[1]
		print ""
		print "%s" %(feature)
		print "%f" %(std)
		print 'What proportion of each race or ethnicity has the specific trait?'
		pprint.pprint(dict(identified_features_info[feature]['proportions']))
		print "How many people in each race had the specific trait?"
		pprint.pprint(dict(identified_features_info[feature]['counts']))
	
	for item in reversed(sorted_by_std):
		print item[0]
	try:
		filename = os.path.join(directory, 'selected.features.p')
		pickle.dump(identified_features_info, open(filename,'wb'))	
	except:
		print "Could not save the selected features info"


	
