from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO

ids_dict = defaultdict(str)
for line in open("Thorax_SampleIDs.txt", "r"):
 if line.startswith("#"):
  continue
 #
 Key,	Pat_ID, Sample_ID = line.rstrip().split()
 ids_dict[Sample_ID] = Key

seq_dict = defaultdict(list)
for line in open("04_THORAX_SABE1172UNHESMSV_ORF_alignment.txt", "r"):
 hit = line.rstrip().split("\t")[-1]
 sampleid, trinity_id, tl, fl, pos, cg, div, pr_seq, og_seq = hit.split("|")
 sid = ids_dict[sampleid] 
 seq_dict[og_seq].append(f"{sid}:{trinity_id}")


outfil = open("07_PredictedORFs_grouped_IDs.txt", "w+")
outfil.write("#TRX_ID\tPopFreq\tAll_IDs\tSeq\tPepSeq\n")
outfil_fa = open("07_PredictedORFs_grouped.fa", "w+")
outfil_fa2 = open("07_PredictedORFs_grouped_cds.fa", "w+")

new_id, count = 0, 0
y = list()

for seq, ids in seq_dict.items():
 new_id += 1
 population = [id.split(":")[1].split("_")[1] for id in ids if "S" not in id.split(":")[0]]
 popfreq = len(set(population)) / 19
 # Look up peptide seq
 sampleid = [i.split(":")[0] for i in ids][0]
 sampleid = [i for i, k in ids_dict.items() if sampleid == k][0]
 seqid =  [i.split(":")[1] for i in ids][0]
 pepsequences = SeqIO.parse(open(f"{sampleid}_trinity.Trinity.fasta.transdecoder.pep"),'fasta')
 pepseq = [str(s.seq) for s in pepsequences if s.id == seqid][0].replace("*", "")
 #
 if popfreq >= 0.75:
  count += 1
 y.append(popfreq)
 #
 outfil.write(f"TRX{new_id}\t{popfreq}\t{','.join(ids)}\t{seq}\t{pepseq}\n")
 outfil_fa.write(f">TRX{new_id}\n{pepseq}")
 outfil_fa2.write(f">TRX{new_id}\n{seq}")

print(f"Total {count} TRX sequences with freq >= 0.75")
print(f"Max pop freq", round(max(y), 2))

fig, ax = plt.subplots()
ax.hist(y, weights=np.ones(len(y))/len(y))
ax.set_title("Predicted ORF Population Frequency in Thorax")
ax.set_xlabel("Population Fraction")
ax.set_ylabel("Predicted ORF fraction")
plt.savefig("07_PredictedORF_PopFreq_Thorax.png")