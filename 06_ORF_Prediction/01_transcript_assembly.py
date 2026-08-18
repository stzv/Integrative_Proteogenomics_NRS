import os
import subprocess
from glob import glob
import sys
from typing import Tuple

trinity_cmd = "trinityrnaseq-v2.13.2/Trinity"


print("Assemble RNA transcripts")

rna_samples = glob("brazil/20220105_thorax_fastq/*_R1_trimmed.fq.gz")


## 

for rna_file in rna_samples:
    #rna_id = "OE1502_T06-12889"
    LEFT = f"/mnt/fedot21/brazil/20220105_thorax_fastq/{rna_id}_R1_trimmed.fq.gz"
    RIGHT = f"/mnt/fedot21/brazil/20220105_thorax_fastq/{rna_id}_R2_trimmed.fq.gz"
    if os.path.isfile(f'{rna_id}_trinity.Trinity.fasta') == False:
        subprocess.call(f"{trinity_cmd} --seqType fq --left {LEFT} --right {RIGHT} --max_memory 5G --output {rna_id}_trinity --full_cleanup --CPU 20", shell = True)


assembled_transcripts = glob(f"*_trinity.Trinity.fasta")
print(f"Assembled {len(assembled_transcripts)} samples")

subprocess.call("email_stepanka.pl Trinity Thorax", shell = True)
print("All done. Have a nice day!")