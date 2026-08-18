
import collections

genmapping_file = "SABE_1172_UNHESMSV_genomemapping/Based_on_reads_newest/05_SABE_1172_UNHESMSV_10x_anchor_match.txt"
#rnaseq_cjs = "SABE_1172_UNHESMSV_RNASeq/06_SABE_1172_UNHESMSV_NRS_CJ_GTF.txt"
rnaseq_cjs = "SABE_1172_UNHESMSV_RNASeq/06_SABE_1172_UNHESMSV_NRS_best_CJ.txt"

GENMAP_dict = collections.defaultdict(dict)

for line in open(genmapping_file, "r"):
  if line.startswith("#"):
    continue
  #
  nrs, PMC, X10x_loc, X10xbc, reads, anchors_matched, all_anchors = line.rstrip("\n").split("\t")
  anchor_chr, anchor_pos = anchors_matched.split(",")[0].split(":")[:2]
  #
  all_anchors = anchors_matched.split(",")
  all_anchors_sorted = sorted(all_anchors, key = lambda x: x.split(":")[1])
  min_a, max_a = all_anchors_sorted[0], all_anchors_sorted[-1]
  anchor_range = ":".join([anchor_chr, min_a.split(":")[1], max_a.split(":")[1]])
  GENMAP_dict[nrs] = {"anchor_chr": anchor_chr.replace("chr", ""), "anchor_pos": int(anchor_pos), "anchor_range": anchor_range}

RNASEQ_dict = collections.defaultdict(dict)
for line in open(rnaseq_cjs, "r"):
  nrs, position= line.rstrip("\n").split("\t")
  if position.startswith("HLA"):
    continue
  rna_chr = position.split("|")[0].split(":")[0]
  rna_start = position.split("|")[0].split(":")[1].split("-")[0]
  rna_end = position.split("|")[0].split(":")[1].split("-")[1]
  rna_pos = (int(rna_start) + int(rna_end) ) / 2
  RNASEQ_dict[nrs] = {"rna_chr": rna_chr, "rna_pos": rna_pos}

outfile = open("09_Mapping_RNASeq_10X_Anchors2.txt", "w")
outfile.write("nrs\tanchor_chr\tanchor_pos\trna_chr\trna_pos\n")
delim = "\t"

count, count_bp = 0, 0

for nrs in RNASEQ_dict.keys():
  ## Get GENMAP
  genmap_nrs = "_".join(nrs.split("_")[:2])
  length = int(nrs.split("_")[-1]) - int(nrs.split("_")[-2]) + 1
  # if not GENMAP, skip
  if not genmap_nrs in GENMAP_dict.keys():
    continue
  #
  genmap_chr = GENMAP_dict[genmap_nrs]["anchor_chr"]
  genmap_pos = GENMAP_dict[genmap_nrs]["anchor_pos"]
  genmap_range = GENMAP_dict[genmap_nrs]["anchor_range"]
  #
  RNA_chr = RNASEQ_dict[nrs]["rna_chr"]
  RNA_pos = RNASEQ_dict[nrs]["rna_pos"]
  #gene = RNASEQ_dict[nrs]["gene"].replace("exon:gene_id ", "").split(":")[0]
  #gene_pos = ":".join(RNASEQ_dict[nrs]["gene"].split(":")[2:-1])
  #
  if genmap_chr == RNA_chr and abs(genmap_pos - RNA_pos) <= 1000:
    outfile.write(f"{delim.join([nrs, str(length), genmap_range, genmap_chr, str(genmap_pos), RNA_chr, str(RNA_pos)])}\n")
    count += 1
    count_bp += length

print("RNASeq vs 10X vs Anchoring", count, round(count_bp/1000000,2), "Mbp")