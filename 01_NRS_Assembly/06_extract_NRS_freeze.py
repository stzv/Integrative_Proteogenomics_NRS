#!/usr/bin/env python

infile = "SABE_1172_UNHESMSV_merged/SABE1172_unmapped.fa"
outfile = "SABE1172_UNHESMSV_NRS_dark_freeze_final.fa"
freq = "05_frequency_results.txt"

###
import regex
import gzip
from Bio import SeqIO
###

print("Loading contigs")
need = list() #Dictionary with all needed contigs (unmapped and partially mapped)

for line in open(freq, "r"):
    ctg = line.split(";")[0]
    need.append(ctg)

print(f"  To be processed: {len(need):,} contigs (with desired coverage)")
##
print("Processing binary map")
partial = dict()
placing = dict()
binary = list()

counter = 0

for line in gzip.open("04_binary_map.txt.gz", "r"):
    line = line.strip(b"\n")
    ctg, info, bitstr = line.split(b"\t")
    ctg = ctg.decode("utf-8")
    info = info.decode("utf-8")
    bitstr = bitstr.decode("utf-8")
    #
    binary.append(ctg)
    # Filter out contigs that do not meet coverage condition
    if ctg not in need: 
        continue
    # Check the binary string
    if (regex.search("^0{200,}", bitstr) and regex.search("0{200,}$", bitstr) and not regex.search("^0+$", bitstr)): # Check if the string starts with at least 200 '0' and has more than 200 '0' in the middle and doesn't contain all '0'
        partial[ctg] = "1_" + str(len(bitstr))# Length of all the string is saved
    elif (regex.search("^0{200,}", bitstr)): # If the string only starts with at least 200 '0'
        partial[ctg] = "1_" + str(len(regex.search("^0{200,}", bitstr)[0])) # Save only length of '0'
        placing[ctg] = info
    elif (regex.search("0{200,}$", bitstr)): # If at least 200 '0 are in the middle
        partial[ctg] = str(regex.search("0{200,}$", bitstr).start()) + "_"+ str(regex.search("0{200,}$", bitstr).start() + len(regex.search("0{200,}$", bitstr)[0]))
        placing[ctg] = info
    else:
        print("ERROR", ctg, bitstr)
        exit

print(f"  Processed partial: {len(partial):,} contigs")

##
print("Processing fasta sequences")

seqio = SeqIO.index(infile, "fasta")
seqio_out = open(outfile, "w+")

ctg2id = dict()
max_len, tot_len = 0, 0
sequences = list()

for record in seqio:
    ctg = seqio[record].id
    if ctg not in need:
        continue
    ctg_str = seqio[record].seq
    s, e = 1, len(ctg_str)
    # if max_len < e:
        # max_len = e
    if ctg in partial:
        s, e = partial[ctg].split("_")
        ctg_str = ctg_str[(int(s)-1):(int(e))]
    tot_len += len(ctg_str)
    ctg2id[ctg] = f"{ctg}_{s}_{e}"
    if ctg in placing:
        seqio_out.write(f">{ctg2id[ctg]} {placing[ctg]}\n{ctg_str}\n")
    else:
        seqio_out.write(f">{ctg2id[ctg]}\n{ctg_str}\n")

# print(f"Printing sequences to outfile {outfile}")
# SeqIO.write(sequences, outfile, "fasta")

print(f"Total:\t{tot_len/1000000:,.2f} Mbps")
print(f"Total:\t{len(ctg2id):,} contigs")


print("\nAll done, have a nice day.")
