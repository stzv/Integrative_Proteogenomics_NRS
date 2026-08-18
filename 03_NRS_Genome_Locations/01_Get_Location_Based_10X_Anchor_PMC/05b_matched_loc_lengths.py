import collections


print("Load NRS lengths")
NRS_len = collections.defaultdict(dict)
for line in open("SABE_1172_UNHESMSV_genotyping/SABE_1172_UNHESMSV_lens_darkfreeze_final.txt", "r"):
      nrs, len = line.rstrip("\n").split("\t")
      nrs = "_".join(nrs.replace(">","").split("_")[:2])
      NRS_len[nrs] = int(len)


##
locs_infile = open("SABE_1172_UNHESMSV_genomemapping/Based_on_reads_newest/05_SABE_1172_UNHESMSV_10x_anchor_match.txt", "r")
total_len = 0

for line in locs_infile:
	if line.startswith("#"):
		continue
	nrs = line.split("\t")[0]
	total_len += NRS_len[nrs]

print("Total Mbps", round(total_len/1000000, 2))