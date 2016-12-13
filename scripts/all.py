import os
import sys

def main():

    """Checks if we have the right number of command line arguments
       and reads them in"""
    if len(sys.argv) < 1:
        print "you must call program as: python ./all.py <path> <cancertype>"
        sys.exit(1)

    filepath = sys.argv[1] # data/BRCA/BRCA
    cancertype = sys.argv[2]

    clinical_file = filepath + "/" + cancertype + ".clinical/" + cancertype + ".clin.merged.picked.txt"
    pval_file = filepath + "_pvals.tsv"

    os.system('python rna.py ' + filepath + ".data.txt")
    os.system('python rna.py ' + clinical_file)
    os.system('python glm.py ' + filepath + " " + cancertype)
    os.system('python fdr.py ' + pval_file)
    os.system('python posthoc.py ' + filepath + " " + cancertype)

if __name__ == "__main__":
    main()
