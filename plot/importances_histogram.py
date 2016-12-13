import matplotlib.pyplot as plt
import numpy as np


#data
# np.random.seed(42)
# data = np.random.rand(5)
class SNP:
	def __init__(self, label):
		self.label = label 
		temp = label.split('_')
		self.gene = label[0]
		self.chr = label[1]
		self.start_pos = label[2]

def read_mini(mini):
	scores = []
	features = []
	featuress = mini.strip().split('\n')
	for feature in featuress:
		splitted = feature.split(' ')
		score = splitted[0]
		name = splitted[1]
		scores.append(float(score))
		features.append(name)
	return (scores,features)

gbm = '''0.002129 pk_chr2_173239482
0.001911 CALD1_chr7_134928898
0.001854 DENND2A_chr7_140602061
0.001847 FAM46C_chr1_117623598
0.001839 U2SURP_chr3_143054978
0.001819 SSTR2_chr17_73170158
0.001807 GDAP1L1_chr20_44264469
0.001739 MED13_chr17_62029991
0.001697 RIMS1_chr6_72248051
0.001674 SUFU_chr10_102550107'''

kich = '''0.036756 SMCO2_chr12_27502000
0.034590 KIAA0430_chr16_15602137
0.033064 DDX24_chr14_94057863
0.030877 DSCAM_chr21_40708597
0.028588 NUDT9_chr4_87458016
0.019116 CTH_chr1_70432189
0.018458 HMHA1_chr19_1073207
0.017718 DCST2_chr1_155030128
0.017555 ZFR2_chr19_3827543
0.017084 RP11428P162_chr10_77731043'''

kirp = '''0.002858 KIAA1109_chr4_122343498
0.002456 TMEM2_chr9_71685763
0.002401 STON2_chr14_81396111
0.002309 TDRD1_chr10_114218472
0.002063 SMPD4_chr2_130173551
0.001976 BPIFB2_chr20_33020377
0.001955 LRP1_chr12_57195738
0.001930 ASS1P4_chrX_3367678
0.001922 LYZL6_chr17_35934840
0.001896 FBXW2_chr9_120764775'''

lgg = '''0.018252 TAS1R2_chr1_18840326
0.017309 PSMD12_chr17_67345801
0.007168 RP11347H156_chr11_50246614
0.007081 CDH12_chr5_21975348
0.006366 AHR_chr7_17339515
0.006103 bP2171C214_chr21_9027104
0.004748 SLC25A51P1_chr6_65788826
0.004567 MYH8_chr17_10396380
0.004270 AP1S1_chr7_101156711
0.004150 ZBTB20_chr3_114351344'''

luad = '''0.002273 BPIFC_chr22_32432544
0.002269 NOVA1_chr14_26480143
0.001993 FAM90A1_chr12_8224782
0.001791 OR2M1P_chr1_248122080
0.001757 CHML_chr1_241634335
0.001753 IGSF1_chrX_131277088
0.001700 LCE3B_chr1_152613929
0.001620 BMPR1APS2_chr11_121362295
0.001514 ACTL6B_chr7_100643326
0.001384 NCAPD3_chr11_134204904
'''

lusc = '''0.003395 FLNB_chr3_58109308
0.002828 KBTBD7_chr13_41192282
0.002823 ADH1C_chr4_99340600
0.002699 KCNRG_chr13_50015532
0.002579 ZNF181_chr19_34739152
0.002576 AP4M1_chr7_100106731
0.002354 LAMC1_chr1_183135088
0.002318 ARMCX2_chrX_101656894
0.002273 HGF_chr7_81705488
0.002100 DZANK1_chr20_18398537'''

mini_list = [gbm, kich, kirp, lgg, luad, lusc]

luad_r ='''0.003378 LRAT
0.002974 PLD6
0.002545 INPP5F
0.002541 OAS3
0.002526 SLC45A1
0.002399 MUT
0.002363 C16orf88
0.002287 METRNL
0.002238 NFYA
0.002198 C3orf65
0.002183 TBC1D25
0.002168 FOXQ1
0.002134 DFFA
0.002022 CREB3
0.001944 RND1
0.001943 C17orf63
0.001931 ATXN1
0.001928 PDCD7
0.001922 APOBEC3B
0.001826 SLC16A8
0.001818 CNNM3
0.001809 PRSS8
0.001801 CDSN
0.001785 KIAA0562
0.001753 RRP7B
0.001733 C3orf62
0.001731 ATRNL1
'''
gbm_r = '''0.002659 PPAP2A
0.002607 HMGCL
0.002476 GSTM1
0.002361 KIAA1804
0.002246 DECR2
0.002072 MOXD1
0.002060 GSDMC
0.002045 BET1
0.002015 ADIPOR2
0.001851 ACAA2
0.001713 IFLTD1
0.001676 SPDYE3
0.001661 LECT1
0.001661 EID2B
0.001633 SHFM1
0.001626 COL8A1
0.001624 C6orf182
0.001549 HS1BP3
0.001527 FOXD4L1
0.001514 TRPC3
0.001514 DSEL
0.001505 PHF19
0.001504 C17orf70
0.001489 LOC255167
0.001488 CASD1
0.001475 C6orf221'''



rna_list = [luad_r, gbm_r]
rna_list = []
for mini in rna_list:
	scores, features = read_mini(mini)
	print scores 
	print features

	fig = plt.figure()
	ax = plt.subplot(111)
	width=0.3
	bins = map(lambda x: x-width/2,range(1,len(scores)+1))
	ax.bar(bins,scores,width=width)
	ax.set_xticks(map(lambda x: x, range(1,len(scores)+1)))
	ax.set_xticklabels(features,rotation=45, rotation_mode="anchor", ha="right")
	fig.tight_layout()
	fig.savefig('importance_histogram')
	plt.show()

# UTS2
# 0.255659
# What proportion of each race or ethnicity has the specific trait?
lala = {'asian': 0.8548387096774194,
 'black or african american': 0.5396825396825397,
 'white': 0.3499420625724218}
# How many people in each race had the specific trait?
# {'american indian or alaska native': 1.0,
#  'asian': 53.0,
#  'black or african american': 102.0,
#  'white': 302.0}

scores=lala.values()
features=lala.keys()

fig = plt.figure()
ax = plt.subplot(111)
width=0.3
bins = map(lambda x: x-width/2,range(1,len(scores)+1))
ax.bar(bins,scores,width=width)
ax.set_xticks(map(lambda x: x, range(1,len(scores)+1)))
ax.set_xticklabels(features,ha="center")
fig.tight_layout()
fig.savefig('importance_histogram')
plt.show()
