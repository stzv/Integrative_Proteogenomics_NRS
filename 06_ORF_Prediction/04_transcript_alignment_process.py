import gzip
from glob import glob
from Bio import SeqIO


filein = glob("*_ORF_SABE1172UNHESMSV_alignment.sam.gz")
divergence_treshold = 0.02 # IM 98%
completeness_treshold = 0.95

print("Join alignments")

unique_nrs = set()
unique_orf = set()

NRS_dict = dict()
for f in filein:  
  rna_id = f.replace("_ORF_SABE1172UNHESMSV_alignment.sam.gz", "")
  for line in gzip.open(f, "rt"):
    line = line.strip("\n")
    # Skip header
    if line.startswith("@"):
      continue
    # Extract param
    transcript, flag, nrs, pos, mapq, cigar = line.split("\t")[:6]
    seq = line.split("\t")[9]
    trans_len = len(seq)
    predicted_sequences = SeqIO.parse(open(f"{rna_id}_trinity.Trinity.fasta.transdecoder.cds"),'fasta')
    ogseq = [str(s.seq) for s in predicted_sequences if s.id == transcript][0]
    # Skip mapping quality < 20
    if int(mapq) < 20:
      continue
    # Skip suplementary alignment
    if int(flag) & 2048: 
      continue
    # Skip if divergence more than 2% (IM 98%)
    divergence = 0
    divergence = [com for com in line.split("\t")[12:] if "de:" in com][0]
    div = float(divergence.split(":")[-1])
    if div > divergence_treshold:
      continue
    # Add to dict
    if nrs not in NRS_dict:
      NRS_dict[nrs] = list()
    NRS_dict[nrs].append(f"{rna_id}|{transcript}|tl:{trans_len}|fl:{flag}|pos:{pos}|cg:{cigar}|{divergence.replace('de:f','div')}|align_seq{seq}|og_seq:{ogseq}")
    unique_nrs.add(nrs)
    unique_orf.add(transcript)

print(f" Total {len(unique_nrs)} unique NRS within ORFs")
print(f" Total {len(unique_orf)} unique ORFs")

print("Print to file")
outfile = open("04_THORAX_SABE1172UNHESMSV_ORF_alignment.txt", "w+")

for nrs in NRS_dict.keys():
  length = int(nrs.split("_")[-1]) - int(nrs.split("_")[-2]) + 1
  #hits = ",".join(NRS_dict[nrs])
  for hit in NRS_dict[nrs]:
    outfile.write(f"{nrs}\t{length}\t{hit}\n")

outfile.close()

####


print("All done.")