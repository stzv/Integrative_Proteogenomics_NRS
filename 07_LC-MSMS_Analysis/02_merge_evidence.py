from collections import defaultdict
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

print("Processing")

TRXids = defaultdict(dict)
MAIN_dict = defaultdict(dict)
x = dict()
orf_sequences = defaultdict(str)

for line in open("00_Predicted_ORFs_grouped_IDs.txt", "r"):
 if line.startswith("#"):
  continue
 #
 TRX_ID, PopFreq, All_IDs, SEQ, PepSeq = line.rstrip().split("\t")
 ids = [id.split(":")[1] for id in All_IDs.split(",")]
 TRXids[TRX_ID] = All_IDs.split(",")
 x[TRX_ID] = float(PopFreq)
 orf_sequences[TRX_ID] = PepSeq

SupORF_files = glob("01_supported_ORF_*.txt")
supported_trx = set()

for fil in SupORF_files:
 sample_id = fil.split("_")[3].replace(".txt", "")
 for line in open(fil, "r"):
  seq, peptides = line.rstrip().split("\t")
  # Find TRX ID
  trx = [trx_id for trx_id, ids in TRXids.items() for i in ids if seq == i.split(":")[1]][0]
  if trx not in MAIN_dict.keys():
   MAIN_dict[trx] = {}
  # Prepare dictionary
  if sample_id not in MAIN_dict[trx].keys():
   MAIN_dict[trx][sample_id] = set()
  # Add the peptide evidence per sample
  for p in peptides.split(","):
   MAIN_dict[trx][sample_id].add(p)

orf_supported_peptides = [orf for orf, evidence in MAIN_dict.items() if evidence]
print(" ORF supported", len(orf_supported_peptides))

print("Into outfiles")

outfile = open("02_Merged_Peptide_Evidence.txt", "w+")
outfile.write("#ORF\tPopFreq\tSupFreq\tCountSupportingSamples\tCountSupportingPeptidesUnique\tAllEvidence\tORFSeq\n")

outfile2 = open("02_Predicted_NotSupported.txt", "w+")
outfile2.write("#ORF\tPopFreq\n")

outfilefa = open("02_List_All_Supported_ORFs.fa", "w+")

y = dict()

for orf, evidence in MAIN_dict.items():
 popfreq = x[orf]
 supporting_samples = len(evidence.keys())
 y[orf] = supporting_samples/19 #Only 19 samples were used for Peptide run
 supporting_peptides = set()
 all_evidence = list()
 for sample, peptides in evidence.items():
  all_evidence.append(f"{sample}:{','.join(peptides)}")
  for p in peptides:
   supporting_peptides.add(p)
 outfile.write(f"{orf}\t{round(popfreq, 4)}\t{round(supporting_samples/19, 4)}\t{supporting_samples}\t{len(supporting_peptides)}\t{'|'.join(all_evidence)}\t{orf_sequences[orf]}\n")
 outfilefa.write(f">{orf}\n{orf_sequences[orf]}\n")

for trx in TRXids.keys():
 if trx not in MAIN_dict.keys():
  outfile2.write(f"{trx}\t{round(x[trx], 4)}\n")
  y[trx] = 0

print("Into graphs")

fig, ax = plt.subplots()
ax.hist(y.values(), weights=np.ones(len(y))/len(y))
ax.set_title("Predicted ORF Supported Frequency in Thorax (histogram)")
ax.set_xlabel("Population Support Fraction")
ax.set_ylabel("Predicted ORF fraction")
plt.savefig("02_PredictedORF_SupportedFrequency_Thorax.png")

df = pd.DataFrame([x, y])
df = df.transpose().rename(columns = {0: "PopFreq", 1: "SupFreq"}).sort_values(by = "PopFreq").dropna()
ms = df.groupby(['PopFreq','SupFreq']).size().reset_index().rename(columns = {0:'count'})

fig, ax = plt.subplots()
ax.scatter(ms["PopFreq"], ms["SupFreq"], s = 1, c = "blue", marker = "x")
ax.scatter(ms["PopFreq"], ms["SupFreq"], s = ms["count"]/2, alpha = 0.5)
ax.set_ylim(top = 1)
ax.set_xlim(right = 1)
ax.set_title("Predicted ORF Supported Frequency vs Popul Freq")
ax.set_xlabel("Population Freq Fraction")
ax.set_ylabel("Supported Freq Fraction")
plt.grid(visible = 1)
plt.savefig("02_PopFreq_vs_SuppFreq_Thorax.png")