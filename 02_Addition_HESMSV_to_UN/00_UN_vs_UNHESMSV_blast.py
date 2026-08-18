import subprocess


fastafile = "SABE_1172_UNHESMSV_genotyping/SABE_1172_UNHESMSV_NRS_dark_freeze_final.fa"
fastafile2 = "SABE_1172_UN_complete/ref_1172/Sao_Paulo_dark_freeze2.fa"

# command = f"makeblastdb -in {fastafile} -dbtype nucl"
# subprocess.run(command, shell = True)

outformat = '"6 qseqid sseqid pident length qstart qend sstart send evalue"'
command = f"blastn -task megablast -perc_identity 95 -outfmt {outformat} -query {fastafile2} -db {fastafile} -out 00_UN_vs_UNHESMSV.m6"
subprocess.run(command, shell = True)

