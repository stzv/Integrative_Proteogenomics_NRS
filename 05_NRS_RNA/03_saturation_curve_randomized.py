import math
from os import access
import random


print("Load in coverage data")
coverage_file = "02_SABE_1172_UNHESMSV_NRS_RNA_reads_count.txt"


coverage_data = dict()
for line in open(coverage_file, "r"):
    if line.startswith("#"):
        accessids = line.strip("\n").split("\t")[1:]
    if line.startswith("k141"):
        ctg = line.strip("\n").split("\t")[0]
        counts = line.strip("\n").split("\t")[1:]
        coverage_data[ctg] = list()
        for c in counts:
            coverage_data[ctg].append(c.split(":")[0]) # Keep only RNAID, no need for reads count

print("Process curve")
sample_count = len(accessids)
random.shuffle(accessids) # randomize the id list
covered_nrs, curve_data, idslist = set(), list(), dict()
cumulative_coverage = 0
NRS_count = 572414

counter = 0
runnum = 3
rna_coverage = open(f"03_SABE_1172_UNHESMSV_RNAcoverage_{sample_count}.txt", "w+")
rna_coverage_addition = open(f"03_SABE_1172_UNHESMSV_RNAcoverageaddition_{sample_count}_run{runnum}.txt", "w+")
previous_cov = 0

for id in accessids:
    subset = [k for k, v in coverage_data.items() if id in v] # Subset NRS covered by RNAID
    #if subset: # if list not empty
    covered_nrs.update(subset) # Add them to the covered NRS list
    curve_data.append(len(covered_nrs)*100/NRS_count)
    coverage = len(covered_nrs)*100/NRS_count
    idslist[id] = coverage
    rna_coverage.write(f"{id}\t{len(subset)*100/NRS_count}\n")
    # Save NRAID's addition to coverage
    addition = (coverage) - previous_cov
    rna_coverage_addition.write(f"{id}\t{addition}\n")
    previous_cov = coverage # update previous
    

print(f"Total coverage {max(curve_data):.0f}%")

with open(f"03_SABE_1172_UNHESMSV_saturation_curve_data_{sample_count}_byRNA_run{runnum}.txt", "w+") as outfile:
    for k in sorted(idslist, key=idslist.get, reverse=True):
        outfile.write(f"{k}\t{idslist[k]}\n")


print("Create graphs")
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.use('Agg')

fig, ax = plt.subplots()

ax.plot(range(1, len(curve_data)+1, 1),curve_data)
ax.set_ylabel("% of NRS Covered")
ax.set_xlabel("RNASeq samples")

fig.tight_layout()
plt.savefig(f"03_SABE_1172_UNHESMSV_RNASeq_saturation_curve_{sample_count}_randomized_byRNA_run{runnum}.png")

import subprocess
subprocess.call("email_stepanka.pl saturation curve", shell = True)
print("All done, have a nice day!")