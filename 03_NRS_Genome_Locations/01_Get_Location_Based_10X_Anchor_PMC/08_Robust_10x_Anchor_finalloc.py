
from asyncio import subprocess
import collections
import sys

print("Load NRS lengths")
NRS_len = collections.defaultdict(dict)
for line in open("SABE_1172_UNHESMSV_genotyping/SABE_1172_UNHESMSV_lens_darkfreeze_final.txt", "r"):
      nrs, len = line.rstrip("\n").split("\t")
      nrs = "_".join(nrs.replace(">","").split("_")[:2])
      NRS_len[nrs] = int(len)

print("Load NRS locations")

outfile = open("08_Robust_mapping_10xPlusAnchor.txt", "w")
NRS_dict = collections.defaultdict(dict)
located_len = 0

for line in open("SABE_1172_UNHESMSV_genomemapping/Based_on_reads_newest/05_SABE_1172_UNHESMSV_10x_anchor_match.txt", "r"):
	if line.startswith("#"):
		continue
	#
	nrs, PMC, X10_loc, X10bc, reads, anchors_matched, all_anchors = line.rstrip("\n").split("\t")
	anchor_start = anchors_matched.split(",")[0]
	anchor_end = anchors_matched.split(",")[-1]
	best_loc = anchor_start.split(":")[0].replace("chr", "") + ":" +str(round(int(anchor_start.split(":")[1]) + int(anchor_end.split(":")[1]) / 2))
	#
	NRS_dict[nrs] = best_loc
	located_len += NRS_len[nrs]

#print(" Total", len(NRS_dict.keys()), "NRS with location")
print("  Which is", round(located_len/1000000), "Mbps")

sys.exit()

print("Load in GTF")

GTF_dict = collections.defaultdict(dict)

for line in open("GTF/Homo_sapiens.GRCh38.105_original.gtf", "r"):
      #Skip header
      if line.startswith("#"):
            continue
      # Split GTF entry
      seqname, source, feature, start, end, score, strand, frame, attributes = line.rstrip("\n").split("\t")
      if not feature == "gene":
            continue
      # Add to dictionary
      gene_id = [attr for attr in attributes.split("; ") if attr.startswith("gene_id")][0].replace("gene_id ", "").replace('"', "")
      gene_name = [attr for attr in attributes.split("; ") if attr.startswith("gene_name")]
      if gene_name:
            gene_name = gene_name[0].replace("gene_name ", "").replace('"', "")
      else:
            gene_name = ""
      gene_biotype = [attr for attr in attributes.split("; ") if attr.startswith("gene_biotype")][0].replace("gene_biotype ", "").replace('"', "").replace(';', "")
      #
      GTF_dict[gene_id] = {"gene_id": gene_id,"gene_name": gene_name, "chr": seqname, "start": start, "end": end, "biotype": gene_biotype, "feature": feature}


print("Process annotation")
outfile = open("08_NRS_Annotation_genes.txt", "w+")
gene_list = set()

for gene, info in GTF_dict.items():
	gene_range = range(int(info["start"]), int(info["end"]))
	for nrs, loc in NRS_dict.items():
		nrs_chr, nrs_pos = loc.split(":")
		if info["chr"] == nrs_chr and int(nrs_pos) in gene_range: # if gene and nrs same chromosome and nrs within gene span
			outfile.write(f'{",".join([nrs, loc, info["gene_name"], info["chr"], info["start"], info["end"]])}\n')
			gene_list.add(f'{gene},{info["chr"]}, {info["start"]},{info["end"]}\n')

outfile2 = open("08_Genes_With_NRS.txt", "w+")

for gene in gene_list:
	outfile2.write(gene)
