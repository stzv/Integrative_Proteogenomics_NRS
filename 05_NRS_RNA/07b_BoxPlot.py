import collections
from statistics import mean
import sys


# data = pd.read_csv("SABE_1172_UNHESMSV_RNASeq/07_SABE_1172_UNHESMSV_NRS_RNA_tissue_table.txt", sep = "\t", header = 0, decimal = ".")
# data['BestTissue'] = data.Tissues.str.split(",").str[0]
# data[['Sample','CPM','TissueType']] = data.BestTissue.str.split(":", expand = True)

# data['CPM'] = data['CPM'].astype(float)

# data = data.loc[data['CPM'] >= 10]

# tissues_examine = ['blood','lung','heart','skin','liver','brain','retina']
# sns.set(style="darkgrid")

# for t in tissues_examine:
# 	data2 = data[data["TissueType"] == t]
# 	average = data2["CPM"].mean()
# 	print(t, average, data2["CPM"].min(),  data2["CPM"].max())

Tissue_dict = collections.defaultdict(dict)

for line in open("SABE_1172_UNHESMSV_RNASeq/07_SABE_1172_UNHESMSV_NRS_RNA_tissue_table.txt", "r"):
	if line.startswith("#"):
		continue
	#
	nrs, popfreq, transfreq, Tissues_Expression, Exons = line.rstrip("\n").split("\t")
	# if nrs == "k141_27904025_1_1077":
	# 	print(line)
	# 	sys.exit()
	if float(transfreq) == 0:
		continue
	#
	tissues_expression_list = Tissues_Expression.split(",")
	for t in tissues_expression_list:
		# Skip errorenous tissue entries
		if len(t.split(":")) < 3:
			continue
		#
		rnasample, cpm, tissue = t.split(":")
		if not float(cpm) > 10:
			continue
		if tissue not in Tissue_dict.keys():
			Tissue_dict[tissue]["CPM"] = list()
			Tissue_dict[tissue]["NRS"] = set()
			Tissue_dict[tissue]["Samples"] = set()
		#
		Tissue_dict[tissue]["CPM"].append(float(cpm))
		Tissue_dict[tissue]["NRS"].add(nrs)
		Tissue_dict[tissue]["Samples"].add(rnasample)

outfile = open("07b_Tissue_Mbps_expression.tab", "w+")
outfile.write(f"tissue\ttissue_expression_mbps\ttissue_expression_meanCPM\ttissue_samples\n")

#print(Tissue_dict["lung"]["NRS"])

delim = '\t'
tissue_expr_sizes = collections.defaultdict(dict)

for key, values in Tissue_dict.items():
	tissue_bps, tissue_mbps = 0, 0
	# Calculate NRS size of expression
	for nrs in values["NRS"]:
		length = int(nrs.split("_")[-1]) - int(nrs.split("_")[-2]) + 1
		tissue_bps += length
	tissue_mbps = round(tissue_bps/1000000, 2)
	if tissue_mbps >= 0.1:
		tissue_expr_sizes[key] = tissue_mbps
	# Calculate averate tissue expression per sample
	average_expression = round(mean(values["CPM"]), 2)
	#print(key, average_expression, len(values["Samples"]))
	#
	outfile.write(f"{delim.join([key, str(tissue_mbps), str(average_expression), str(len(values['Samples'])), ','.join(values['Samples'])])}\n")


	
