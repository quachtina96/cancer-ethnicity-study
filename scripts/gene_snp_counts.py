'''
@fileoverview This file extracts the tumor sammple barcodes and corresponding 
hugo symbols for each snp within a given folder containing MAF files respresenting
SNP data from TCGA (The Cancer Genome Atlas) on GDC (Genomic Data Commons).

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
            Hugo_Symbol
            Tumor_Sample_Barcode"""
    column_indices = [0,15]
    snp_data = np.genfromtxt(filename, delimiter='\t', dtype=None, skip_header=1, names=True, usecols=column_indices)
    return snp_data

def getPatientBarcodeHugoSymbolPairs(maf_file_list):
    """Maps the first four terms of the corresponding
    Tumor_Sample_Barcode to the Hugo_Symbol."""

    def writeToFile(snp_data, filename):
        hugo = snp_data['Hugo_Symbol']
        barcodes = snp_data["Tumor_Sample_Barcode"]
        f = open(filename, 'w')
        for i in xrange(len(hugo)):
            hugo_symbol = hugo[i]
            barcode = '-'.join(barcodes[i].split('-')[:4])
            # Use the first four sections of the barcode.
            # barcode_to_uuid_map[barcode] = hugo_symbol
	    f.write('\t'.join([hugo_symbol,barcode]))
	    f.write('\n')
    barcode_hugo_files = []
    for maf in maf_file_list:
        snp_data = readMaf(maf)
        # mapPatientBarcodeToUUID(barcode_to_uuid_map, snp_data)
        # mapBarcodeToSnp(sample_to_snps, snp_data, snpSet)
	
	maf_filename = maf.split('/')[-1]
	print maf_filename
        maf_id = maf_filename.split('.')[:4]
	print maf_id
        filename = 'barcode_hugo_'+ '.'.join(maf_id) + '.txt'
        writeToFile(snp_data, filename)
        barcode_hugo_files.append(filename)
    return barcode_hugo_files

if __name__ == '__main__':
    # parse command-line arguments
    if len(sys.argv) < 1:
        print "you must call program as:  "
        print "   python snp.py <TUMOR_TYPE_DIR>"
        sys.exit(1)

    tumor_type_directory = sys.argv[1]

    maf_list = []
    print "List of MAF files using for SNP analysis..."
    for root, dirs, files in os.walk(tumor_type_directory):
        for file in files:
            if file.endswith("maf.txt"):
                maf_list.append(os.path.join(root, file))
                print os.path.join(root, file)

    print "Processing MAF files..."
    filenames = getPatientBarcodeHugoSymbolPairs(maf_list)
    
    for filename in filenames:
	os.system("mv " + filename + " " + tumor_type_directory + "/" + filename)
    
    print "Files mapping patient barcode to Hugo Symbol can be found in %s" %(tumor_type_directory)
