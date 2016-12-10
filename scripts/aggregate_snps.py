import numpy as np
import csv
import sys
import pickle


def aggregate_snps(filename):
	snp_info = {}
	with open(filename) as f:
		for line in f:
			parts = line.strip().split('\t') # Gene PatientID
			print parts
			if parts[0] in snp_info:
				snp_info[parts[0]].append(parts[1])
			else:
				snp_info[parts[0]] = [parts[1]]
	pickle.dump(snp_info, open(filename + '.dict.p', 'wb'))

def main():

    """Checks if we have the right number of command line arguments
       and reads them in"""
    if len(sys.argv) < 1:
        print "you must call program as: python ./aggregate_snps.py <datafile>"
        sys.exit(1)

    filename = sys.argv[1]
    aggregate_snps(filename)

if __name__ == "__main__":
    main()

