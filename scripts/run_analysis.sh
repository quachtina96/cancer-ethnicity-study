#!/bin/bash

# python randomforest.py ../data/gdc_data/BRAIN/GBM > ../results/randomforest/TCGA.GBM.RandomForest.txt
# python randomforest.py ../data/gdc_data/BRAIN/LGG > ../results/randomforest/TCGA.LGG.RandomForest.txt
python randomforest.py ../data/gdc_data/LUNG/LUAD > ../results/randomforest/TCGA.LUAD.RandomForest.txt
python randomforest.py ../data/gdc_data/LUNG/LUSC > ../results/randomforest/TCGA.LUSC.RandomForest.txt
python randomforest.py ../data/gdc_data/BRCA > ../results/randomforest/TCGA.BRCA.RandomForest.txt
python randomforest.py ../data/gdc_data/KIDNEY/KICH > ../results/randomforest/TCGA.KICH.RandomForest.txt
python randomforest.py ../data/gdc_data/KIDNEY/KIRP > ../results/randomforest/TCGA.KIRP.RandomForest.txt
python randomforest.py ../data/gdc_data/KIDNEY/KIRC > ../results/randomforest/TCGA.KIRC.RandomForest.txt


