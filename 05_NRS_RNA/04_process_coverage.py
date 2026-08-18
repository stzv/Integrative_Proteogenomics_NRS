import subprocess
import numpy as np

coverage_file = "02_SABE_1172_UNHESMSV_NRS_RNA_reads_count.txt"
NRS_fasta = "../01_NRS_Assembly/SABE_1172_UNHESMSV_NRS_dark_freeze_final.fa"
#NRS_fasta = "SABE1172_UNHESMSV_NRS_dark_freeze_final.fa"

print("Load in NRS library")
NRS_dict = dict()
for line in open(NRS_fasta, "r"):
    if line.startswith(">"):
        ctg = line.strip("\n").split(" ")[0].strip(">")
        length = int(line.strip("\n").split(" ")[0].split("_")[3]) - int(line.strip("\n").split("_")[2]) + 1
        NRS_dict[ctg] = length

print(f" {len(NRS_dict)} NRS")

print("Load in coverage data")
CVG_dict = dict()
RNAIDs_dict = dict()
RNAcount = []
for line in open(coverage_file, "r"):
    line = line.strip("\n")
    #
    if line.startswith("#"):
        RNAcount = line.split("\t")[1:]
    if line.startswith("k141"):
        ctg = line.split("\t")[0]
        counts = [c.split(":")[1] for c in line.split("\t")[1:]]
        rnaids = line.split("\t")[1:]
        if counts:
                CVG_dict[ctg] = counts
                RNAIDs_dict[ctg] = rnaids

print(f" Processed {len(RNAcount)} RNA samples")
print(f" {len(CVG_dict.keys())} NRS with coverage")
print(f" {len(CVG_dict.keys())*100/len(NRS_dict):.2f}% of total NRS")

print("Calculate coverage frequency")
frequency_dict = dict()
for nrs in NRS_dict:
    if nrs in CVG_dict.keys():
        hits = CVG_dict[nrs]
        cov = len(hits)
        #
        frequency_dict[nrs] = (cov / len(RNAcount))

# Sort frequency dict descending
frequency_dict = sorted(frequency_dict.items(), key=lambda x:x[1])
frequency_dict = dict(frequency_dict)

delim = "\t"

print("Print out frequency into file")
with open(f"04_SABE_1172_UNHESMSV_RNASeq_freq.txt", "w+") as outfile:
    outfile.write(f"#NRS\t{delim.join(RNAcount)}\n")
    for nrs in frequency_dict.keys():
        outfile.write(f"{nrs}\t{frequency_dict[nrs]}\t{','.join(RNAIDs_dict[nrs])}\n")


print("Print out frequency cut-off for commonly expressed NRS")
cutoff = 0.1 
with open(f"04_SABE_1172_UNHESMSV_RNASeq_freq_cutoff{cutoff*100:.0f}perc.txt", "w+") as outfile:
    for nrs in frequency_dict.keys():
        if frequency_dict[nrs] >= cutoff:
            outfile.write(f"{nrs}\t{frequency_dict[nrs]}\t{','.join(CVG_dict[nrs])}\n")

print("Create graphs")
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
import numpy as np

fig, ax = plt.subplots()

allele_freq_list = np.array([v for v in frequency_dict.values() if v >= 0.1])

ax.hist(allele_freq_list, bins=40, weights=np.ones(len(allele_freq_list))/len(allele_freq_list), color = "lightgreen", edgecolor='darkgreen', linewidth=0.5)
#ax.set_ylim(0, 0.1)
ax.set_xlim(0.1, 1)
ax.set_ylabel("Fraction of contigs")
ax.set_xlabel("RNASeq frequency")

fig.tight_layout()
plt.savefig("04_SABE_1172_UNHESMSV_RNASeq_freq_hist.png")

#subprocess.call("email_stepanka.pl process coverage", shell = True)
print("\nAll done, have a nice day!")