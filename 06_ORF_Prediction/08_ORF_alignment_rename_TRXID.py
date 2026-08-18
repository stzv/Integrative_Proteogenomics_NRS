import collections


TRXids = collections.defaultdict(dict)
for line in open("07_PredictedORFs_grouped_IDs", "r"):
 if line.startswith("#"):
  continue
 #
 TRX_ID,	PopFreq,	All_IDs,	SEQ, PepSeq = line.rstrip().split("\t")
 ids = [id.split(":")[1] for id in All_IDs.split(",")]
 TRXids[TRX_ID] = ids

outfile = open("08_THORAX_SABE1172UNHESMSV_ORF_alignment_TRXIDs.txt", "w+")

for line in open("04_THORAX_SABE1172UNHESMSV_ORF_alignment.txt", "r"):
 nrs, nrs_len, alignment = line.rstrip().split("\t")
 #OE1507_04-1252|TRINITY_DN1010_c0_g1_i1.p1|tl:336|fl:0|pos:154|cg:230M106S|div:0|ATGGC
 orfid = alignment.split("|")[1]
 trxid = [key for key, value in TRXids.items() if orfid in value][0]
 outfile.write(f"{nrs}\t{nrs_len}\t{trxid}\t{alignment}\n")
