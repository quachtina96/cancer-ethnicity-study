'''
@fileoverview This file handles iterative random forest analysis on cancer
data from TCGA (The Cancer Genome Atlas) on GDC (Genomic Data Commons).

@author Tina Quach (quacht@mit.edu)
'''

import numpy as np 
from randomforest import RandomForest
from collections import defaultdict
from feature_importance import FeatureImportances

def iterative_random_forest(data, class_labels, n_iterations, leave_out_proportion):
	'''Given data as a numpy array (whose shape is (n_samples, n_features) 
	excluding the labels) and corresponding feature_labels as a numpy array with
	shape (n_features,), fit the random forest iteratively, removing 20% of the
	least significant features based on importance.'''
	rf = RandomForest()
	fi = FeatureImportances()

	# 0 means the values are valid
	label_mask = np.array([0]*data.shape[1])
	data_mask = np.array([[0]*data.shape[1]]*data.shape[0])
	masked_data = np.ma.masked_array(data, data_mask)

	# Initial fit
	init_fit = rf.fit(data, class_labels)
	fi.set_importances(rf.classifier.feature_importances_)

	classifiers = [init_fit]
	oob_scores =[init_fit.oob_score_]
	for i in xrange(n_iterations):
		# Remove the from the data 20% least important features 
		n_features = np.sum(label_mask)
		n_leave_out = round(n_features*leave_out_proportion)
		
		least_important = fi.get_least_important_features_indices(n_leave_out)

		for index in least_important:
			label_mask[index] = 1
			data_mask[:,index] = np.ones((1,data.shape[0]))

		masked_data = np.ma.masked_array(data, data_mask)
		newest_fit = rf.fit(masked_data, class_labels)
		print 'iteration %d' %(i)
		classifiers.append(newest_fit)
		oob_scores.append(rf.classifier.oob_score_)
		print 'oob_score: %f' %(rf.classifier.oob_score_)
		fi.set_importances(rf.classifier.feature_importances_)

	print 'Classifiers'
	print classifiers
	best_classifier = classifiers[np.argmax(oob_scores)]
	return best_classifier, classifiers

if __name__ == '__main__':
	# TODO(generalize the loading of the SNP matrix)
	data = np.load('LGG_snp_matrix.tsv.npy')
	data = np.array([list(line) for line in data])
	snp_data = data[:,1:-2]
	class_labels = data[:,-2]
	feature_labels = np.load('LGG_snp_matrix.tsv_snps.npy')
	n_iterations = 4
	leave_out_proportion = .20
	optimal_random_forest = iterative_random_forest(data, class_labels, n_iterations, leave_out_proportion)

