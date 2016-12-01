'''
@fileoverview This file handles the single nucleotide variation analysis on
data from TCGA (The Cancer Genome Atlas) on GDC (Genomic Data Commons).

@author Tina Quach (quacht@mit.edu)
'''

import numpy as np 
from collections import defaultdict
import demographic
import os
import subprocess

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
#
# @param{numpy.ndarray} an array of tuples with the following attributes in 
# order.
	# Chromosome
	# Start_Position
	# Tumor_Sample_Barcode
# @param{defaultdict<defaultdict>} sample_to_SNPs - defaultdict mapping samples 
# to SNPs
# @param{set<tuple>} - snpSet 
def mapSampleToSnp(sample_to_Snps, snp_data, snpSet):
	for snp in snp_data:
		snp_loc = (snp["Chromosome"], snp["Start_Position"])
		sample_to_Snps[snp["Tumor_Sample_Barcode"]][snp_loc] += 1
		snpSet.add(snp_loc)
	#return sample_to_Snps

# Filter SNPs and return a dictionary mapping samples to SNPs that meet the 
# threshold.
def filterSampleSnpMap(sample_to_snps, threshold):
	filteredSnpSet = set()
	filteredSampleMap = defaultdict(defaultdict)
	if threshold == 1:
		# Avoid unecessary computation. No filtering will occur.
		return

	for barcode in sample_to_snps:
		filtered = {snp: count for snp, count in sample_to_snps[barcode].iteritems() if count >= threshold}
		print filtered.values()

		# Only include sample if its SNPs meet the threshold
		if filtered:
			filteredSampleMap[barcode] = filtered
			for snp in sample_to_snps[barcode]:
				filteredSnpSet.add(snp)
	return filteredSampleMap

# Given a set of possible SNPs and map of samples to SNPs, create a file of a 
# given format that represents a matrix of binary values representing whether or 
# not a given patient/sample contains the SNP.
def formMatrix(snpSet, sample_to_snps, filename, format):
	# Open file for writing
	pass

# Using BRCA data, in a test data folder outside of the codebase,
# applied the functions above in order to filter the data
def test():
	# Let's assume that I'm given the names of all the files relevant to a specific
	# cohort or TCGA project...as a list
	path = '../../gdc_data/test_BRCA/'

	# Get the working directory and strip off new line character so that we can 
	# switch back to this original directory. Then, change directory to the one
	# containing desired files
	original_working_dir  = subprocess.check_output("pwd", shell=True).rstrip()
	os.chdir(path)

	file_list = ['TCGA.BRCA.somaticsniper.09fe60b3-0dd8-411d-ad26-69a81410fb98.somatic.maf.txt',
	'TCGA.BRCA.varscan.608b713a-e27d-4b68-8728-9b12c500c077.somatic.maf.txt',
	'TCGA.BRCA.muse.8751a889-cb3e-4487-ba6f-ac91651666e7.somatic.maf.txt',
	'TCGA.BRCA.mutect.9408fdf2-013f-4c09-8821-a709af56b9ff.somatic.maf.txt']

	# Given a list of MAF files for a specific cohort, creates a dictionary mapping 
	# each sample to the SNPs found in that sample.
	int_default_dict = lambda: defaultdict(int)
	sample_to_snps = defaultdict(int_default_dict)
	snpSet = set()
	for maf in file_list:
		snp_data = readSimpleMaf(path + maf)
		mapSampleToSnp(sample_to_snps, snp_data, snpSet)

	print len(snpSet)
	tumor_sample_barcodes = sample_to_snps.keys()

	# Filter SNPs and return a dictionary mapping samples to SNPs that meet the 
	# threshold.
	threshold = len(file_list)/4 + 1
	filtered = filterSampleSnpMap(sample_to_snps, threshold)
	print len(filtered)

	os.chdir(original_working_dir.rstrip())




