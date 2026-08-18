import subprocess

samfile = "GRCh38_UNHESMSV_vs_Chm13v2.sam.gz"

#subprocess.run(f"gunzip {samfile} -k", shell = True)
subprocess.run(f"cut -f1 {samfile.replace('.gz', '')} | grep k141 | uniq > ctg_ids.txt", shell = True)

outfile = open("ctg_lengths.txt", "w")
outfile.write("NRS,Length\n")

for line in open("ctg_ids.txt", "r"):
    nrs = line.strip("\n")
    short_nrs = "_".join(nrs.split("_")[:2])
    length = int(nrs.split("_")[-1]) - int(nrs.split("_")[-2]) + 1
    outfile.write(",".join([nrs, short_nrs, str(length)]) + "\n")
