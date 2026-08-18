from glob import glob
import subprocess
import os

fastafile = "SABE_1172_UNHESMSV_genotyping/SABE1172_UNHESMSV_NRS_dark_freeze_final.fa" 
transcripts = glob(f"*_trinity.Trinity.fasta.transdecoder.cds")
#transcripts = ["OE1502_MRC5-PS_trinity.Trinity.fasta.transdecoder.cds"]

for trans_file in transcripts:
 trans_id = trans_file.replace("_trinity.Trinity.fasta.transdecoder.cds", "")
 #
 subprocess.call(f"minimap2 -ax sr -t 24 --sam-hit-only {fastafile} {trans_file} -o {trans_id}_ORF_SABE1172UNHESMSV_alignment.sam; gzip {trans_id}_ORF_SABE1172UNHESMSV_alignment.sam", shell = True)

#
subprocess.call("email_stepanka.pl THORAX ORF NRS alignment", shell = True)