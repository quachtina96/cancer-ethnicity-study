#!/bin/bash

# SNP DATA ANALYSIS
# python randomforest.py ../data/gdc_data/LUNG/LUAD snp 10 > ../results/randomforest/TCGA.LUAD.RandomForest_2.txt
# python randomforest.py ../data/gdc_data/LUNG/LUSC snp 10 > ../results/randomforest/TCGA.LUSC.RandomForest_2.txt
# python randomforest.py ../data/gdc_data/BRCA snp 10 > ../results/randomforest/TCGA.BRCA.RandomForest_2.txt
# python randomforest.py ../data/gdc_data/KIDNEY/KICH snp 10 > ../results/randomforest/TCGA.KICH.RandomForest_2.txt
# python randomforest.py ../data/gdc_data/KIDNEY/KIRP snp 10 > ../results/randomforest/TCGA.KIRP.RandomForest_2.txt
# python randomforest.py ../data/gdc_data/KIDNEY/KIRC snp 10 > ../results/randomforest/TCGA.KIRC.RandomForest_2.txt
# python randomforest.py ../data/gdc_data/BRAIN/GBM snp 10 > ../results/randomforest/TCGA.GBM.RandomForest_2.txt
# python randomforest.py ../data/gdc_data/BRAIN/LGG snp 10 > ../results/randomforest/TCGA.LGG.RandomForest_2.txt

# RNA DATA ANALYSIS
python randomforest.py /home/cancer-ethnicity-study/data/rna_data/BRCA rna 15 > /home/cancer-ethnicity-study/results/randomforest/rna/BRCA.rf.results.txt
python randomforest.py /home/cancer-ethnicity-study/data/rna_data/LUAD rna 15 > /home/cancer-ethnicity-study/results/randomforest/rna/LUAD.rf.results.txt
python randomforest.py /home/cancer-ethnicity-study/data/rna_data/LGG rna 15 > /home/cancer-ethnicity-study/results/randomforest/rna/LGG.rf.results.txt

# Analyze Feature Importances


