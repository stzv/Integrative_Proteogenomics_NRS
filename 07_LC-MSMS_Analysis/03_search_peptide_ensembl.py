from Bio import SeqIO
import gzip
import collections
from fnmatch import fnmatch

import re
import sys

##
print(" Load ensembl database")
ensembl = SeqIO.parse(gzip.open("/home/szverinova/ensembl_hs_pep/Homo_sapiens.GRCh38.pep.all.fa.gz", "rt"),'fasta')

##
print(" Load peptides")
pep_set = set()
peptides =  collections.defaultdict(dict)
for line in open("01_peptides_list.txt", "r"):
   if line.startswith("#"):
      continue
   p = line.rstrip().split("\t")[0]
   freq = int(line.rstrip().split("\t")[1])/19
   peps = line.rstrip().split("\t")[-3]
   pep_set.add(p)
   peptides[p] = [freq, peps]

##
print(" Load spectrums per peptide")
spectrums = collections.defaultdict(dict)
for line in open("01_peptides_list.txt", "r"):
   pep = line.rstrip().split("\t")[0]
   ids  = line.rstrip().split("\t")[-1]
   spectrums[pep] = ids

##
print(" Check for known peptides")
known_peptides = set()
database_count = 0
annotated_peptides = collections.defaultdict(dict)

for en in ensembl:
   database_count += 1
   for p in pep_set:
      # Skip found peptide if it is already known
      if p in known_peptides:
         continue
      # Create regex search for (Iso)Leucin
      p_il = re.sub("[LI]", "(L|I)", p)
      # Match found peptide to ensembl entry
      match = re.search(p_il, str(en.seq))
      if match:
         # Save the original format of found peptide if known in ensembl
         known_peptides.add(p)
         if p not in annotated_peptides.keys():
            annotated_peptides[p] = list()
         annotated_peptides[p].append(str(en.seq))
   # Counter
   if database_count % 10000 == 0:
      print(database_count)

print(" Database contains", database_count, "entries")


print("\n Total peptides checked", len(pep_set))
print(" Peptides found", len(known_peptides))

orf_sequences = collections.defaultdict(str)
for line in open("00_Predicted_ORFs_grouped_IDs.txt", "r"):
   if line.startswith("#"):
      continue
   #
   TRX_ID, PopFreq, All_IDs, SEQ, PepSeq = line.rstrip().split("\t")
   ids = [id.split(":")[1] for id in All_IDs.split(",")]
   for i in ids:
      orf_sequences[i] = PopFreq

## Print found peptides
outpep2 = open("03_Peptides_Found_Ensembl.txt", "w+")
outpep2.write("#Peptide\tEvidenceFrequence\tEvidenceSeq\tAnnotatedPeptides\n")
annotated_peptides_count = set()

for p, annotated in annotated_peptides.items():
   for ap in annotated:
      annotated_peptides_count.add(ap)
   #
   outpep2.write(f"{p}\t{peptides[p][0]}\t{peptides[p][1]}\t{','.join(annotated)}\n")

print(" Found total of", len(annotated_peptides_count), "annotated peptides containing evidence peptides")

# Print not found peptides
TRXids = collections.defaultdict(dict)

for line in open("00_Predicted_ORFs_grouped_IDs.txt", "r"):
 if line.startswith("#"):
  continue
 #
 TRX_ID,	PopFreq,	All_IDs,	SEQ, PepSeq = line.rstrip().split("\t")
 ids = [id.split(":")[1] for id in All_IDs.split(",")]
 TRXids[TRX_ID] = ids

outpep = open("03_Peptides_NotFound_Ensembl.txt", "w+")
outpep.write("#Peptide\tPredictedSeqFreq\tEvidenceFrequence\tTRXIds\tOrfs\tMSSpectrums\n")
unknown_count = 0

for p in pep_set:
   if p not in known_peptides:
      unknown_count += 1
      seqs = peptides[p][1]
      predfreq = list()
      trxid = set()
      for s in seqs.split(","):
         predfreq.append(orf_sequences[s])
         for key, value in TRXids.items():
            if s in value:
               trxid.add(key)
      outpep.write(f"{p}\t{','.join(predfreq)}\t{peptides[p][0]}\t{','.join(trxid)}\t{seqs}\t{spectrums[p]}\n")

print(" Peptides not found in database", unknown_count)

import subprocess
subprocess.call("email_stepanka.pl Peptide search ensembl", shell = True)
