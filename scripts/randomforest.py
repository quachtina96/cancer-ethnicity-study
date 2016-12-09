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

def readSNPMatrix(matrix_file):
	# Read the Patient_Barcode, first 4 SNP data, the Race, and the Ethnicity
	# fields.
	with open(matrix_file, 'r') as f:
		labels = f.readline().strip().split('\t')

	data = np.genfromtxt(matrix_file, delimiter='\t', dtype=None,names=True)
	return (labels, data)

def random_forest(matrix_path):
	# matrix_path = '../../gdc_data/BRAIN/GBM/GBM.tsv'
	with open(matrix_path, 'r') as f:
		labels = f.readline().strip().split('\t')

	col_indices = range(1,len(labels)-2)
	print "Getting Patient Barcodes..."
	patient_barcodes = np.genfromtxt(matrix_path, delimiter='\t',usecols=[0], dtype=None,names=True)
	print "Reading in SNP Data..."
	data = np.array(np.genfromtxt(matrix_path, delimiter='\t',usecols=col_indices, dtype=None,names=True))
	race = np.genfromtxt(matrix_path, delimiter='\t',usecols=[len(labels)-2], dtype=None,names=True)   
	ethnicity = np.genfromtxt(matrix_path, delimiter='\t',usecols=[len(labels)-1], dtype=None,names=True)   
	random_forest = RandomForestClassifier(n_estimators=15, oob_score=True)
	fitted_random_forest = random_forest.fit(data, race.reshape(-1, 1))
	return fitted_random_forest

# if __name__ == '__main__':
#     # parse command-line arguments
 #    if len(sys.argv) < 1:
 #        print "you must call program as:  "
 #        print "   python randomforest.py <matrix_dir>"
 #        sys.exit(1)
 #    directory = sys.argv[1]

 #    matrix_list = []
 #    print "List of Matrix TSV files using for Random Forest analysis..."
 #    for root, dirs, files in os.walk(directory):
 #        for file in files:
 #            if file.endswith(".tsv"):
 #            	print file
 #                matrix_list.append(os.path.join(root, file))
 #                print os.path.join(root, file)
    
 #    matrix_path = matrix_list[0]
 #    matrix_filename = matrix.split('/')[-1]
    
 #    with open(matrix_path, 'r') as f:
	# 	labels = f.readline().strip().split('\t')

	# col_indices = range(1,len(labels)-2)
	# patient_barcodes = np.genfromtxt(matrix_file, delimiter='\t',usecols=[0], dtype=None,names=True)
	# data = np.genfromtxt(matrix_file, delimiter='\t',usecols=col_indices, dtype=None,names=True)   
	# race = np.genfromtxt(matrix_file, delimiter='\t',usecols=[len(labels-2)], dtype=None,names=True)   
	# ethnicity = np.genfromtxt(matrix_file, delimiter='\t',usecols=[len(labels-1)], dtype=None,names=True)   

 #    random_forest = RandomForestClassifier(n_estimators=15, oob_score=True)
	# fitted_random_forest = random_forest.fit(data, race)
