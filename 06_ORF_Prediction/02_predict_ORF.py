import subprocess
from glob import glob
import os

import sys

TransAssembly_files = glob("*_trinity.Trinity.fasta")

for file in TransAssembly_files:
 rnaid = file.replace("_trinity.Trinity.fasta", "")
 # Run Transdecoder on the sequence
 if os.path.isfile(f"{rnaid}_trinity.Trinity.fasta.transdecoder.cds") == False:
  subprocess.call(f"TransDecoder.LongOrfs -t {file}", shell = True)
  subprocess.call(f"TransDecoder.Predict -t {file}", shell = True)
 # Remove intermediate files
 #subprocess.call(f"rm -r 04_{rnaid}.fasta.transdecoder_dir*")

####
print("All done, have a good day!")
subprocess.call("email_stepanka.pl ORF prediction", shell = True)