'''
@fileoverview This file handles the single nucleotide variation analysis on
data from TCGA (The Cancer Genome Atlas) on GDC (Genomic Data Commons).

@author Tina Quach (quacht@mit.edu)
'''

import numpy as np 
from collections import defaultdict
import os
import subprocess
import sys

def readSNPMatrix(matrix_file):
	# Read the Patient_Barcode, first 4 SNP data, the Race, and the Ethnicity
	# fields.
	with open(matrix_file, 'r') as f:
		labels = f.readline().strip().split('\t')
	data = np.genfromtxt(matrix_file, delimiter='\t', dtype=None,names=True)
	return (labels, data)


# TODO: implement pca (can build off of 047 pset)