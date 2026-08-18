from glob import glob
from Bio import SeqIO
import subprocess

sampleids_file = open("Thorax_SampleIDs.txt", "r")
sampleids_dict = dict()

for line in sampleids_file:
 if line.startswith("#"):
  continue
 #
 patient_id, sampleid, fileid = line.rstrip().split()
 sampleids_dict[fileid] = patient_id

# filelist = glob("*_trinity.Trinity.fasta")
# for f in filelist:
#  if f.replace("_trinity.Trinity.fasta", "") not in sampleids_dict.keys():
#   print(f.replace("_trinity.Trinity.fasta", ""))


orf_file = open("04_THORAX_SABE1172UNHESMSV_ORF_alignment.txt", "r")
sample_dict = dict()

for line in orf_file:
 nrs, nrslen, transcript = line.rstrip().split("\t")
 rna_id, trans_id, trans_len, flag, pos, cigar, div, trans_seq = transcript.split("|")
 #
 sample_dict.setdefault(rna_id, list())
 sample_dict[rna_id].append(trans_id)

counter = 0

for rna, transcripts in sample_dict.items():
 counter += 1
 print(counter)
 if rna not in sampleids_dict.keys():
  continue
 #
 patid = sampleids_dict[rna]
 #
 outfasta = open(f"Predicted_ORF_inclNRS/{patid}_ORF.fasta", "w+")
 pep_file = SeqIO.parse(open(f"{rna}_trinity.Trinity.fasta.transdecoder.pep"),'fasta')
 #
 for fa in pep_file:
  name, seq = fa.id, str(fa.seq)
  if name in sample_dict[rna]:
   outfasta.write(f">{name}\n{seq.replace('*', '')}\n")
 outfasta.close()

subprocess.call("tar -zcvf Predicted_ORF_inclNRS.tar.gz Predicted_ORF_inclNRS/", shell = True)