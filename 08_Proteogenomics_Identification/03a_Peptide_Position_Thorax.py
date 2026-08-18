from Bio import SeqIO
import collections
import re
import math

import sys

####
print("Load in orf & pep sequences")

ORF_dict = collections.defaultdict(dict)
for record in SeqIO.parse("Thorax/07_PredictedORFs_grouped_cds.fa", "fasta"):
    ORF_dict[str(record.id).split()[0].strip(">")] = str(record.seq)

PEP_dict = collections.defaultdict(dict)
for record in SeqIO.parse("Thorax/07_PredictedORFs_grouped.fa", "fasta"):
#     for nrs in NRS_dict.keys():
#         if nrs in record.id:
    PEP_dict[str(record.id).split()[0].strip(">")] = str(record.seq)

####
print("Load in alignment information")
NRS_dict = collections.defaultdict(list)
for line in open("Thorax/00_ORF_Alignment.txt", "r"):
    #k141_31560909_1_693	693	TRX1_726|0.1928374655647383|fl:16|pos:1|cg:66S140M520S|div:0
    nrs, nrs_len, alignment = line.rstrip().split("\t")
    orfid, al, fl, pos, cg, div = alignment.split("|")
    orfid = orfid.split("_")[0]
    # Get NRS loc coord on peptide seq
    cg_split = re.findall("\d+\w", cg.replace("cg:", ""))
    if int(fl.replace("fl:", "")) == 16:
        cg_split = cg_split[::-1]
    nrs_overlap = ""
    for cg in cg_split:
        if cg[-1] == "S":
            nrs_overlap = nrs_overlap + "0" * int(cg[:-1])
        elif cg[-1] in ["M", "I"]:
            nrs_overlap = nrs_overlap + "1" * int(cg[:-1])
        elif cg[-1] == "D": # Deletion -> NRS extra information that is missing in reference => Skip
            continue
    nrs_start = nrs_overlap.find("1")
    nrs_end = nrs_overlap.rfind("1")
    nrs_loc = f"{nrs_start}:{nrs_end}"
    nrs_start_pep = int(nrs_start/3)
    nrs_end_pep = int(nrs_end/3)
    ## Cannot get completely precise nt to pep position due to possible framshift => int rounding allows for 1 aa error
    nrs_loc_pep = f"{nrs_start_pep}:{nrs_end_pep}"
    # # Add values
    value = {"nrs": nrs, "pep_location": nrs_loc_pep}
    if value not in NRS_dict[orfid]:
        NRS_dict[orfid].append(value)


####
print("Check Peptide positions on NRS")
outfile = open("03a_Peptide_Position_on_ORF.tsv", "w+")
outfile.write("#Peptide\tNRS\tWithinNRS\tPep_Pos\tORF_Pos\tnew_ORF\tnew_AA\n")

outfile2 = open("03a_NovelPeptides_NRS_Thorax.txt", "w+")

for line in open("03_Peptides_NotFound_Ensembl_Thorax.txt", "r"):
    if line.startswith("#"):
        outfile2.write(line)
        continue
    ##
    Peptide = line.split("\t")[0]
    ORF = line.split("\t")[4]
    ##
    print(ORF)
    #print(PEP_dict[ORF])
    AASeq = PEP_dict[ORF]
    ORFSeq = ORF_dict[ORF]
    NRS_PEP_loc = NRS_dict[ORF]
    ##
    pep_start = AASeq.index(Peptide)
    pep_end = AASeq.index(Peptide) + len(Peptide)
    Pep_Pos = str(pep_start + 1) + ":" + str(pep_end)
    ## Pep overlap with NRS    
    within_nrs = "False"
    for pep_loc in NRS_PEP_loc:
        nrs_pep_start, nrs_pep_end = pep_loc["pep_location"].split(":")
        if int(nrs_pep_start) <= pep_start and pep_end <= int(nrs_pep_end):
            within_nrs = "True"
    ##
    pep_start_orf = pep_start*3
    pep_end_orf = pep_end*3
    ORF_Pos = str(pep_start_orf + 1) + ":" + str(pep_end_orf)
    ##
    new_AA = AASeq[:pep_start] + AASeq[pep_start:pep_end].lower() + AASeq[pep_end:]
    new_ORF = ORFSeq[:pep_start_orf] + ORFSeq[pep_start_orf:pep_end_orf].lower() + ORFSeq[pep_end_orf:]
    ##
    if within_nrs == "True":
        outline = "\t".join([Peptide, ORF, Pep_Pos, ORF_Pos, new_ORF, new_AA])
        outfile.write(f"{outline}\n")
        outfile2.write(line)

