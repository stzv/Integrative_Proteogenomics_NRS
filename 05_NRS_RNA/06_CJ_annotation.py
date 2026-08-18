import gzip
import os
from typing import DefaultDict


####
print("Load NRS with sufficient coverage")
covered_nrs = set()
for line in open("SABE_1172_UNHESMSV_RNASeq/04_SABE_1172_UNHESMSV_RNASeq_freq.txt", "r"):
    nrs = line.split()[0]
    covered_nrs.add(nrs)

####
print("Load CJs")
cj_dict = dict()


counter = 0
cj_dict = dict()
rna_dict = dict()

for line in open("05_RNASeq_chimeric_junctions_merged_binned.txt", "r"):
  nrs, bins, rna = line.strip("\n").split("\t")
  if nrs not in covered_nrs:
    continue
  # Add nrs to dictionary
  cj_dict.setdefault(nrs, {})
  # Add count for each CJ bin
  bins = bins.split(",")
  for b in bins:
    if b.split("|")[0] not in cj_dict[nrs]:
      cj_dict[nrs].setdefault(b.split("|")[0], int(b.split("|")[1]))
    else:
      cj_dict[nrs][b.split("|")[0]] += int(b.split("|")[1])


print("Count final CJ")
NRS_CJ_file = open("06_SABE_1172_UNHESMSV_NRS_best_CJ.txt", "w+")
final_cj_dict = dict()

for nrs in cj_dict.keys():
    maximum = max([count for bucket, count in cj_dict[nrs].items()])
    if maximum < 5:
        continue
    final_cj_dict[nrs] = [f"{bucket}|{count}" for bucket, count in cj_dict[nrs].items() if int(count) == maximum]
    NRS_CJ_file.write(f"{nrs}\t{','.join(final_cj_dict[nrs])}\n")

####
print("Load in GTF Annotation file")
gene_dict = dict()
human_chr = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"]

for line in gzip.open("/mnt/fedot21/brazil/20210622_STAR/Homo_sapiens.GRCh38.104.gtf.gz", "r"):
 # Skip header
    if line.startswith(b"#"):
        continue
    # Split line into values
    chrname, source, feature, start, end, score, strand, frame, attributes = line.decode().strip("\n").split("\t")
    # Skip non-human chromosomes (decoys, ....)
    if chrname not in human_chr:
        continue
    # Create new entry for gene if not already existing in dictionary
    #gene_id = attributes.split(";")[0].strip(" gene_id ")
    exon_id = attributes.split(";")[0].strip(" exon_id ")
    if chrname not in gene_dict.keys():
        gene_dict[chrname] = list()
    # Save to dictionary
    #if feature == "exon":
    entry = {"Chr":chrname, "Start": int(start), "End": int(end), "Str":strand, "Feat":feature, "attributes":attributes}
    if entry not in gene_dict[chrname]:
        gene_dict[chrname].append(entry)

####
print("Check CJs against GTF")
print(f" Total {len(final_cj_dict.keys())} NRS with CJ")
gtf_outfile = open("06_SABE_1172_UNHESMSV_NRS_CJ_GTF_full.txt", "w+")
counter = 0

for nrs in final_cj_dict:
    counter += 1
    if counter % 100000 == 0:
        print(counter)
    for cj in final_cj_dict[nrs]:
        # Skip CJ if not on human chromosome or MT
        if cj.split(":")[0] not in human_chr:
            continue
        # Extract GTF entries on the chromosome as the CJ
        gtf_list = gene_dict[cj.split(":")[0]]
        # Keep GTF entries where CJ falls within it's range
        entry_list = [g for g in gtf_list if (int(cj.split("|")[0].split(":")[1].split("-")[0]) >= g["Start"] and int(cj.split("|")[0].split(":")[1].split("-")[1]) <= g["End"])]
        # Save the GTF entries for each CJ of each NRS
        for g in entry_list:
            gtf_outfile.write(f'{nrs}\t{cj}\t{g["Feat"]}:{g["Chr"]}:{g["Start"]}:{g["End"]}:{g["Str"]}:{g["attributes"]}\n')
    #

print("\nAll done, have a nice day!")