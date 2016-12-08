import csv

rna_patients = set()
snp_patients = set()
with open('../BRCA_RNA_patients.csv', 'rb') as csvfile:
	reader = csv.reader(csvfile)
	for row in reader:
		for item in row:
			rna_patients.add('-'.join(item.split('-')[:4]))

with open('barcodeSet.txt', 'rb') as snpfile:
	for line in snpfile:
		snp_patients.add('-'.join(line.strip().split('-')[:4]))

# print rna_patients ^ snp_patients
count = 0
for item in rna_patients:
	if item in snp_patients:
		print item
		count += 1

print len(rna_patients), 'rna' # need to filter out normal samples
print len(snp_patients), 'snp'
