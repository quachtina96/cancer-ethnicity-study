
# BRCA: TRABD 18286, BAX 1375, CTU2 4405, PWP2 14369
# GBM: HOTAIR, FOXE1
# LGG: ATP2A1, C7orf31, BTBD16, ARMC10, CFHR1, LAIR2, C12orf5
# LUAD:
# LUSC:
def build_thing(cancer, num, gene):
	print 'sudo python scripts/posthoc.py data/rna_data/' + cancer + "/" + cancer + " " + cancer
	print 'sudo python scripts/heatmap.py data/rna_data/' + cancer + "/" + cancer + " " + cancer
	print 'sudo python scripts/boxplot.py data/rna_data/' + cancer + "/" + cancer + " " + cancer + " " + str(num) + " " + gene

# build_thing('BRCA', 18286, 'TRABD')
# build_thing('BRCA', 1375, 'BAX')
# build_thing('BRCA', 4405, 'CTU2')
# build_thing('BRCA', 14369, 'PWP2')
build_thing('GBM', 1, 'der')