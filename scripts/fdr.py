from statsmodels.sandbox.stats.multicomp import fdrcorrection0, multipletests
import sys
import math
# Input filename of genes + pvalues
def fdr_correct(filename):
	genes = []
	pvals = []
	with open(filename) as f:
		for line in f:
			parts = line.strip().split('\t')
			gene = parts[0]
			pval = float(parts[1])

			if math.isnan(pval):
				print gene + ' has a p-value of ' + parts[1]
			else:
				genes.append(gene)
				pvals.append(pval)

	corrected = multipletests(pvals, method='fdr_bh')
	new_pvals = corrected[1]

	outfile = open('.'.join(filename.split('.')[:-1]) + '.corrected.tsv', 'w')

	accepted = 0
	for i, val in enumerate(new_pvals):
		if val < 0.05:
			accepted += 1
			outfile.write(str(genes[i]) + "\t" + str(val) + "\n")

	print str(accepted) + '/' + str(len(new_pvals)) + ' genes accepted.'


def main():

    """Checks if we have the right number of command line arguments
       and reads them in"""
    if len(sys.argv) < 1:
        print "you must call program as: python ./fdr.py <pvalue_file>"
        sys.exit(1)

    filename = sys.argv[1]
    fdr_correct(filename)

if __name__ == "__main__":
    main()
#fdr_correct('../data/BRCA/BRCA_pvals.tsv')
   