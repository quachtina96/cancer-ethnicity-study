'''
@fileoverview This file handles the single nucleotide variation analysis on
data from TCGA (The Cancer Genome Atlas) on GDC (Genomic Data Commons).

@author Tina Quach (quacht@mit.edu)
'''

import numpy as np 
from collections import defaultdict
import demographic

# Read desired field from the Single Nucleotide Variation MAF file 
# (as specified by GDC)
#
# @param{string} filename - name of the maf file
# @return{numpy.ndarray} an array of tuples with the following attributes in 
# order.
	# Hugo_Symbol
	# Entrez_Gene_Id
	# Chromosome
	# Start_Position
	# Variant_Classification
	# Reference_Allele
	# Tumor_Seq_Allele1
	# Tumor_Seq_Allele2
	# Tumor_Sample_Barcode
	# IMPACT
def readMaf(filename):
	column_indices = [0,1,4,5,8,10,11,12,15,92]
	snp_data = np.genfromtxt(filename, delimiter='\t', dtype=None, names=True, skip_header=1, usecols=column_indices)
	return snp_data

# Read desired field from the Single Nucleotide Variation MAF file 
# (as specified by GDC)
#
# @param{string} filename - name of the maf file
# @return{numpy.ndarray} an array of tuples with the following attributes in 
# order.
	# Chromosome
	# Start_Position
	# Tumor_Sample_Barcode
def readSimpleMaf(filename):
	column_indices = [4,5,15]
	snp_data = np.genfromtxt(filename, delimiter='\t', dtype=None, names=True, skip_header=1, usecols=column_indices)
	return snp_data


# Map tumor Sample Barcodes to the SNPs that were found from in the tumor 
# encoded by the barcode.
# @param{dict} sample_to_SNPs - defaultdict mapping samples to SNPs
def mapSampleToSNP(sample_to_SNPs, snp_data):

	sample_to_SNPs = defaultdict(list)
	for snp in snp_data:
		sample_to_SNPs[snp["Tumor_Sample_Barcode"]].append(snp)
	return sample_to_SNPs

# Let's assume that I'm given the names of all the files relevant to a specific
# cohort or TCGA project...as a list

