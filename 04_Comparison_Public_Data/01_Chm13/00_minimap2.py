import subprocess

reference = "GRCh38/GRCh38_SABE_1172_UNHESMSV_extended_ref.fa"
fasta = "chm13v2.0.fa"

command = f"minimap2 -ax sr -t 24 {fasta} {reference} | gzip -c > GRCh38_UNHESMSV_vs_Chm13v2.sam.gz"
subprocess.run(command, shell=True)
