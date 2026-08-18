import os
from Bio import SeqIO
import collections
import numpy as np
from matplotlib import pyplot as plt

print("Load in reads coverage")


unhesmsv_readcount = "SABE_UNHESMSV_readcount_per_nrs.txt"
NRS_UNHESMSV_coverage = collections.defaultdict(int)
if os.path.isfile(unhesmsv_readcount) == False:
    for line in open("01_NRS_Assembly/SABE_1172_UNHESMSV_NRS_dark_freeze_final_coverage.txt", "r"):
        if line.startswith("Contig"):
            continue
        ctg = line.split("\t")[0]
        counts = line.strip("\n").split("\t")[1:]
        for count in counts:
            NRS_UNHESMSV_coverage[ctg] = max(NRS_UNHESMSV_coverage[ctg], int(count))
    with open(unhesmsv_readcount, "w+") as outfile:
        for nrs, readcount in NRS_UNHESMSV_coverage.items():
            outfile.write(f"{nrs}\t{readcount}\n")

un_readcount = "SABE_UN_readcount_per_nrs.txt"
NRS_UN_coverage = collections.defaultdict(int)
if os.path.isfile(un_readcount) == False:
    for line in open("Sao_Paulo_dark_freeze_coverage2.txt", "r"):
        if line.startswith("Contig"):
            continue
        ctg = line.split("\t")[0]
        counts = line.strip("\n").split("\t")[1:]
        for count in counts:
            NRS_UN_coverage[ctg] = max(NRS_UN_coverage[ctg], int(count))
    with open(un_readcount, "w+") as outfile:
        for nrs, readcount in NRS_UN_coverage.items():
            outfile.write(f"{nrs}\t{readcount}\n")


x_ax = [0, 200, 400, 600, 800, 1000, 2000, 5000, 10000, 50000, 2000000]

def LoadCoverage(fil):
    count_dict = dict()
    ctg_set = set()
    ## Get the counts for each read coverage bucket
    for line in open(fil, "r"):
        ctg, count = line.strip("\n").split("\t")
        ctg_set.add(ctg)
        for i in range(len(x_ax)-1):
            #
            if x_ax[i] not in count_dict.keys():
                count_dict[x_ax[i]] = 0
            #
            if int(count) >= x_ax[i] and int(count) < x_ax[i+1]:
                count_dict[x_ax[i]] += 1
    ## Sort the keys
    sort_key = list(count_dict.keys())
    sort_key.sort()
    count_list_sorted = {i: count_dict[i] for i in sort_key}
    ## Recalculate coverage count to frequency
    total_count = len(ctg_set)
    freq_list_sorted = {i: round(count_list_sorted[i]/total_count, 2) for i in count_list_sorted}
    ##
    return list(freq_list_sorted.values())


NRS_UNHESMSV_coverage = LoadCoverage(unhesmsv_readcount)
NRS_UN_coverage = LoadCoverage(un_readcount)

labels = ["0", "200", "400", "600", "800", "1'000", "2'000", "5'000", "10'000", "50'000+"]
index = np.arange(0, len(labels))


####
print("Create the graph")

fig, ax = plt.subplots(figsize=((5.9, 5)))

ax.bar(index + 0.05, NRS_UNHESMSV_coverage, align = 'edge', width = 0.4, color = "b", edgecolor = "b", alpha = 0.2, linewidth = 0.5, label = "UNHESMSV")
ax.bar(index + 0.5, NRS_UN_coverage, align = 'edge', width = 0.4, color = "r", edgecolor = "r", alpha = 0.2, linewidth = 0.5, label = "UN")

ax.set_xticks(np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]))
ax.set_xticklabels(labels, rotation = 20)
ax.set_xlabel('Reads Coverage', size = 11, labelpad = 0.5)
ax.set_ylabel('Fraction of contigs', size = 11)

ax.tick_params(axis="x", labelsize=9)
ax.tick_params(axis="y", labelsize=9) 

ax.set_xlim([0, 9])

#ax.legend()

fig.tight_layout(pad=1)
plt.savefig("01b_UN_vs_UNHESMSV_coverage_v2.png")

## Sanity check
print(np.mean(NRS_UNHESMSV_coverage))