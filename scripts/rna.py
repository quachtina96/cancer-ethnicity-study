import numpy as np
import csv
import sys

# Each column is a patient ID
# First column is a gene
# Up to line 31 are unknown genes

def save_rna_matrix(filename):
	patients = None
	index = 0
	end_range = -1
	genes = []
	num_skip = None

	with open(filename) as f:
		for line in f:
			if index == 0:
				patients = line.split("\t")[1:]
				end_range = len(patients) + 1

			if index > 1:
				parts = line.split("\t")
				gene_id = parts[0]
				if gene_id[0] != "?":
					if num_skip == None:
						num_skip = index
					genes.append(gene_id.split('|')[0])
			index += 1

	data = np.loadtxt(filename, delimiter='\t', skiprows=num_skip, usecols=range(1,end_range))

	out = open(filename + '.patients.csv', 'wb')
	wr = csv.writer(out, quoting=csv.QUOTE_ALL)
	wr.writerow(patients)

	gene_out = open(filename + '.genes.csv', 'wb')
	gr = csv.writer(gene_out, quoting=csv.QUOTE_ALL)
	gr.writerow(genes)

	# Transpose matrix (row is patient, column is gene)
	#data = data.transpose()

	assert len(genes) == data.shape[0]
	assert len(patients) == data.shape[1]

	np.save(filename + '.matrix', data)


def main():

    """Checks if we have the right number of command line arguments
       and reads them in"""
    if len(sys.argv) < 1:
        print "you must call program as: python ./rna.py <datafile>"
        sys.exit(1)

    filename = sys.argv[1]
    save_rna_matrix(filename)

if __name__ == "__main__":
    main()

