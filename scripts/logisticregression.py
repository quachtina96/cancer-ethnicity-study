'''
@fileoverview This file handles the logistic regression on cancer
data from TCGA (The Cancer Genome Atlas) on GDC (Genomic Data Commons).

@author Tina Quach (quacht@mit.edu)
'''

import numpy as np 
from collections import defaultdict
import os
import subprocess
import sys
from sklearn.linear_model import LogisticRegression

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
	
#TODO(quacht): code the logistic regression analysis