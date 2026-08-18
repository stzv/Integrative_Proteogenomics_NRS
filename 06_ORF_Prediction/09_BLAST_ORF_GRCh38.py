import subprocess

outformat = '"6 qseqid sseqid pident length qstart qend sstart send evalue"'
ref = "GRCh38/GRCh38_SABE_1172_UNHESMSV_extended_ref.fa"

#blastn = "ncbi-blast-2.13.0+/bin/blastn"

command = f"blastn -task megablast -perc_identity 95 -outfmt {outformat} -query 07_PredictedORFs_grouped_cds.fa -db {ref} -out 09_Thorax_ORF_GRh38_extended_IM95.txt"
subprocess.call(command, shell = True)
