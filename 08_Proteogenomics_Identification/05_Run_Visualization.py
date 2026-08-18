import subprocess

for line in open("04_NovelPeptides_Merged.tsv", "r"):
	if line.startswith("#"): continue
	peptide = line.split("\t")[0]
	print(peptide)
	subprocess.call(f"Rscript 05_Visualization.R {peptide}", shell = True)