#!/bin/bash


# RNA DATA ANALYSIS
# python randomforest.py ../data/rna_data/BRCA rna 50 1 > ../results/randomforest/rna/_BRCA.rna.rf.results.txt
# python randomforest.py ../data/rna_data/LGG rna 50 1 > ../results/randomforest/rna/_LGG.rna.rf.results.txt
# python randomforest.py ../data/rna_data/LUAD rna 50 1 > ../results/randomforest/rna/_LUAD.rna.rf.results.txt
# python randomforest.py ../data/rna_data/LUSC rna 50 1 > ../results/randomforest/rna/_LUSC.rna.rf.results.txt
# python randomforest.py ../data/rna_data/GBM rna 50 1 > ../results/randomforest/rna/_GBM.rna.rf.results.txt


# SNP DATA ANALYSIS
python randomforest.py ../data/gdc_data/LUNG/LUAD snp 50 1 > ../results/randomforest/_LUAD.snp.rf.results.txt
python randomforest.py ../data/gdc_data/LUNG/LUSC snp 50 1 > ../results/randomforest/_LUSC.snp.rf.results.txt
python randomforest.py ../data/gdc_data/BRCA snp 50 1 > ../results/randomforest/_BRCA.snp.rf.results.txt
python randomforest.py ../data/gdc_data/BRAIN/GBM snp 50 1 > ../results/randomforest/_GBM.snp.rf.results.txt
python randomforest.py ../data/gdc_data/BRAIN/LGG snp 50 1 > ../results/randomforest/_LGG.snp.rf.results.txt

