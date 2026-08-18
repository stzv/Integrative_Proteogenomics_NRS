#!/usr/bin/env python
from Bio import SeqIO

import sys

print("Loading dark freeze sequences")
dark_freeze = set()

with open("SABE1172_UNHESMSV_NRS_dark_freeze_final.fa") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        ctgid = record.id.split("_")[0] + "_" + record.id.split("_")[1]
        dark_freeze.add(ctgid)

print(f" Total sequences {len(dark_freeze)}")

samples = set()
table = dict()
counter = 0

print("Loading in counts dictionary")
for line in open('05_frequency_results.txt', "r"):
    ctg = line.split(";")[0]
    ## Remove after corrected version of 05_run -> no headers are printed there
    if "Contig" in ctg:
        continue
    ##
    if ctg not in dark_freeze:
        continue
    ##
    counter += 1
    counts = line.replace("\n", "").split(";")[1:]
    table[ctg] = {}
    for c in counts:
        if "total" in c: 
            continue
        sampleid, count = c.split(":")
        table[ctg][sampleid] = int(count)
        samples.add(sampleid)

print(f"Total {len(samples)} samples")
print(f"Total {counter} contigs")

print("Iterating through samples")
for contig, vals in table.items():
    for s in samples:
        if s not in vals:
            table[contig][s] = 0

print("\tPrinting the table") #To make sure that every contig has the same order of samples
outfile = open("SABE1172_UNHESMSV_NRS_dark_freeze_final_coverage.txt", "w+")
delim = "\t"
outfile.write(f"Contig\t{delim.join(samples)}\n")


for contig, values in table.items():
    table[contig] = {k: table[contig][k] for k in samples}

for contig, values in table.items():
    counts = delim.join(str(v) for c, v in values.items())
    outfile.write(f'{contig}\t{counts}\n')

del(table, dark_freeze, samples)
print("All done, have a nice day")
