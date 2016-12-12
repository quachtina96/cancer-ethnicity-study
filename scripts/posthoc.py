from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison
import sys
import pandas as pd
from glm import load_data

def tukey(filename, cancer):
	data, patients, genes, master_df, clinical, indices_to_delete = load_data(filename, cancer)
	items = []
	outfile = open(filename + '.posthoc.txt', 'w')
	with open(filename + '_pvals.corrected.tsv', 'r') as f:
		for line in f:
			parts = line.strip().split('\t')
			ind = int(parts[0])
			gene = parts[1]
			pval = float(parts[2])

			items.append((ind, gene, pval))

	sorted_items = sorted(items,key=lambda x:(x[2],x[0]))
	print 'Items sorted by p-values'
	for tup in sorted_items[:20]:
		ind = tup[0]
		expression_data = [item for a, item in enumerate(list(data[ind])) if a not in indices_to_delete]
		df = master_df.copy()
		df['expression'] = pd.Series(expression_data, index=master_df.index)
		assert len(df['expression']) == len(df['gender']) == len(df['stage']) == len(df['race']) == len(df['age'])

		mc = MultiComparison(df['expression'], df['race'])
		result = mc.tukeyhsd()
		outfile.write(str(tup[1]) + '\n')
		outfile.write(str(result) + '\n')
		outfile.write(str(mc.groupsunique) + '\n' + '\n')
	print 'Post-hoc output available at ' + filename + '.posthoc.txt'

def main():

    """Checks if we have the right number of command line arguments
       and reads them in"""
    if len(sys.argv) < 1:
        print "you must call program as: python ./posthoc.py <datafile> <cancertype>"
        sys.exit(1)

    filename = sys.argv[1]
    cancer_type = sys.argv[2]
    tukey(filename, cancer_type)

if __name__ == "__main__":
    main()

#tukey('../data/BRCA/BRCA', 'BRCA')

