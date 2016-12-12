#!/bin/bash

# python randomforest.py ../data/gdc_data/BRAIN/GBM > ../results/randomforest/TCGA.GBM.RandomForest.txt
# python randomforest.py ../data/gdc_data/BRAIN/LGG > ../results/randomforest/TCGA.LGG.RandomForest.txt
python randomforest.py ../data/gdc_data/LUNG/LUAD snp 10 > ../results/randomforest/TCGA.LUAD.RandomForest_2.txt
python randomforest.py ../data/gdc_data/LUNG/LUSC snp 10 > ../results/randomforest/TCGA.LUSC.RandomForest_2.txt
python randomforest.py ../data/gdc_data/BRCA snp 10 > ../results/randomforest/TCGA.BRCA.RandomForest_2.txt
python randomforest.py ../data/gdc_data/KIDNEY/KICH snp 10 > ../results/randomforest/TCGA.KICH.RandomForest_2.txt
python randomforest.py ../data/gdc_data/KIDNEY/KIRP snp 10 > ../results/randomforest/TCGA.KIRP.RandomForest_2.txt
python randomforest.py ../data/gdc_data/KIDNEY/KIRC snp 10 > ../results/randomforest/TCGA.KIRC.RandomForest_2.txt
python randomforest.py ../data/gdc_data/BRAIN/GBM snp 10 > ../results/randomforest/TCGA.GBM.RandomForest_2.txt
python randomforest.py ../data/gdc_data/BRAIN/LGG snp 10 > ../results/randomforest/TCGA.LGG.RandomForest_2.txt

# TODO(quacht): random forest usage has changed

# SNP DATA ANALYSIS

# RNA DATA ANALYSIS

# Analyze Feature Importances


