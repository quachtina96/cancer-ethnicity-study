#!/bin/bash 
# Make SNP Matrices

echo 'BRCA'
python snp.py /home/cancer-ethnicity-study/data/gdc_data/BRCA
echo 'LGG'
python snp.py /home/cancer-ethnicity-study/data/gdc_data/BRAIN/LGG
#echo 'GBM'
#python snp.py /home/cancer-ethnicity-study/data/gdc_data/BRAIN/GBM
echo 'LUAD'
python snp.py /home/cancer-ethnicity-study/data/gdc_data/LUNG/LUAD
#echo 'LUSC'
#python snp.py /home/cancer-ethnicity-study/data/gdc_data/LUNG/LUSC

