import statsmodels.sandbox.stats.multicomp.fdrcorrection0 as fdrcorrection

# Input filename of genes + pvalues
def fdr_correct(filename):
	genes = []
	pvals = []
	with open(filename) as f:
		for line in f:
			parts = line.strip().split('\t')
			genes.append(parts[0])
			pvals.append(parts[1])