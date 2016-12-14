
# BRCA: TRABD 18286, BAX 1375, CTU2 4405, PWP2 14369
def build_thing(cancer, num, gene):
	print 'sudo python scripts/posthoc.py data/rna_data/' + cancer + "/" + cancer + " " + cancer
	print 'sudo python scripts/heatmap.py data/rna_data/' + cancer + "/" + cancer + " " + cancer
	print 'sudo python scripts/boxplot.py data/rna_data/' + cancer + "/" + cancer + " " + cancer + " " + str(num) + " " + gene

build_thing('BRCA', 1375, 'BAX')
build_thing('BRCA', 4405, 'CTU2')
build_thing('BRCA', 14369, 'PWP2')