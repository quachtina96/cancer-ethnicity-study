import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
import seaborn as sns
import csv
from glm import load_data
import sys


def boxplot(filename, cancer, snp_index, gene):
	data, patients, genes, master_df, clinical, indices_to_delete = load_data(filename, cancer)
	expression_data = [item for a, item in enumerate(list(data[snp_index])) if a not in indices_to_delete]
	df = master_df.copy()
	df['expression'] = pd.Series(expression_data, index=master_df.index)
	ax = sns.boxplot(x="race", y="expression", data=df)
	ax.set_title("Gene Expression Comparison for " + gene + " for " + cancer + " patient data")
	sns.plt.savefig(filename + "_boxplot.png")

def main():

    """Checks if we have the right number of command line arguments
       and reads them in"""
    if len(sys.argv) < 1:
        print "you must call program as: python ./glm.py <datafile> <cancertype>"
        sys.exit(1)

    filename = sys.argv[1]
    cancer_type = sys.argv[2]
    snp_index = int(sys.argv[3])
    gene = sys.argv[4]
    boxplot(filename, cancer_type, snp_index, gene)

if __name__ == "__main__":
    main()


# filename = '../data/BRCA/BRCA'
# cancer = 'BRCA'
# snp_index = 18286
# gene = 'TRABD'



