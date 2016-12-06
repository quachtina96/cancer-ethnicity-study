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
import sys

def readMaf(filename):
    """Read desired field from the Single Nucleotide Variation MAF file 
    (as specified by GDC)

    Retrieves rows pertaining to the given keys from the Table instance
    represented by big_table.  Silly things may happen if
    other_silly_variable is not None.

    Args:
        filename: string name of the MAF file

    Returns:
        An array of tuples with the following attributes in order:
            Chromosome
            Start_Position
            Tumor_Sample_Barcode
            Tumor_Sample_UUID"""
    column_indices = [4,5,15,32]
    snp_data = np.genfromtxt(filename, delimiter='\t', dtype=None, names=True, skip_header=1, usecols=column_indices)
    return snp_data

def mapToSnp(dictionary, snp_data, opt_snpSet):
    """Updates the given dictionary by mapping tumor sample UUID to the SNPs
    that were found  in the tumor encoded by the barcode.

    Args:
        dictionary: default dict of default dicts to update with the identifier 
            and SNPs.
        snp_data: a 1-D numpy array of tuples with the following attributes in 
            order:
                Chromosome
                Start_Position
                Tumor_Sample_Barcode
                Tumor_Sample_UUID"""
    for snp in snp_data:
        snp_loc = (snp[0], snp[1]) # (Chromosome, Start_Position)
        uuid = snp[3]
        dictionary[uuid][snp_loc] += 1
        if type(opt_snpSet) == set:
            opt_snpSet.add(snp_loc)


def filterSampleSnpMap(sample_to_snps, threshold):
    """Filter SNPs by the number of times they appear. Return a dictionary 
    mapping tumor sample identifiers to SNPs that meet the given threshold.

    Args:
        dictionary: default dictionary to update with the identifier and SNPs.
        threshold: integer representing the minimum number of times the SNP must
            have been identified in that specific tumor sample
    
    Returns:
        A tuple containing
            1) dictionary mapping tumor sample identifier to the SNPs that are 
            associated with that specific sample and that have been indentified in 
            the given tumor sample at least the threshold number of times    
            2) set of SNPs contained in the filtered dictionary"""
    filteredSnpSet = set()
    filteredSampleMap = defaultdict(defaultdict)
    if threshold == 1:
        # Avoid unecessary computation. No filtering will occur.
        return

    for barcode in sample_to_snps:
        filtered = {snp: count for snp, count in sample_to_snps[barcode].iteritems() if count >= threshold}

        # Only include sample if its SNPs meet the threshold
        if filtered:
            filteredSampleMap[barcode] = filtered
            for snp in sample_to_snps[barcode]:
                filteredSnpSet.add(snp)
    return (filteredSampleMap, filteredSnpSet)

def formMatrix(snpSet, sample_to_snps, empty_file):
    """Given a set of possible SNPs and map of samples to SNPs, create a file of a 
    given format that represents a matrix of binary values representing whether or 
    not a given patient/sample contains the SNP.

    Args:
        snpSet: set of all possible SNPs
        sample_to_snps: map of tumor samples to the SNPs found in them.
        empty_file: name of empty file to write matrix into.
    Returns:
        Numpy array of numpy arrays represnting the matrix.
    """
    snps = list(snpSet)
    # Ensure that the file is empty
    output_file = open(empty_file, 'a')
    if os.stat(empty_file).st_size != 0:
        raise ValueError("Filename provided does not refer to empty file.")
 
    for identifier in sample_to_snps: 
        print identifier
        row = [identifier]
        for snp in snps:
            if snp in sample_to_snps[identifier]:
                row.append(1)
                print "hi"
            else:
                row.append(0)
                print "bye"
        # row += [1 if snp in sample_to_snps[identifier] else 0 for snp in snps]
        demographicInfo = demographic.getDemographicFromTumorUuid(identifier)
        race = demographicInfo.values()[0]
        ethnicity = demographicInfo.values()[1]
        row.append(race)
        row.append(ethnicity)
        output_file.write('\t'.join([str(feature) for feature in row]))
        output_file.write('\n')
    output_file.close()

def processMAF(maf_file_list, matrix_filename):
    """Using BRCA data, in a test data folder outside of the codebase,
    applied the functions above in order to filter the data.

    Args:
        maf_file_list: list of the strings representing the names of the MAF files 
            containing the SNP data to process.
        matrix_filename: string name of the file to write the matrix (result
            of the data processing)"""

    # Given a list of MAF files for a specific cohort, creates a dictionary mapping 
    # each sample to the SNPs found in that sample.
    int_default_dict = lambda: defaultdict(int)
    sample_to_snps = defaultdict(int_default_dict)
    
    snpSet = set()
    print "Reading in SNP MAF files..."
    for maf in maf_file_list:
        snp_data = readMaf(maf)
        mapToSnp(sample_to_snps, snp_data, snpSet)

    # Filter SNPs and return a dictionary mapping samples to SNPs that meet the 
    # threshold.
    print "Filtering SNPs..."
    threshold = len(maf_file_list)/4 + 1
    filteredSampleMap, filteredSnpSet = filterSampleSnpMap(sample_to_snps, threshold)
    
    print "Original Sample Count %d" %(len(sample_to_snps.keys()))
    print "Original SNP Count %d" %(len(snpSet))

    print "Filtered Sample Count %d" %(len(filteredSampleMap.keys()))
    print "Filtered SNP Count %d" %(len(filteredSnpSet))

    print "Forming Matrix..."
    matrix = formMatrix(filteredSnpSet, filteredSampleMap, matrix_filename)

    
    return matrix

if __name__ == '__main__':
    # parse command-line arguments
    if len(sys.argv) < 1:
        print "you must call program as:  "
        print "   python ps1-dotplot.py <TUMOR_TYPE_DIR>"
        sys.exit(1)
    tumor_type_directory = sys.argv[1]

    maf_list = []
    print "List of MAF files using for SNP analysis..."
    for root, dirs, files in os.walk(tumor_type_directory):
        for file in files:
            if file.endswith("maf.txt"):
                maf_list.append(os.path.join(root, file))
                print os.path.join(root, file)
    
    matrix_filename = tumor_type_directory.split('/')[-1] + '.tsv'
    print "Processing MAF files..."
    processMAF(maf_list, matrix_filename)
    
    print "Matrix can be found at %s" %(matrix_filename)




