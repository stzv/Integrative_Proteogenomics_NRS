#!/usr/bin/env python

### ENVIRONMENT VARIABLES ####
min_disc = 0.05
min_unmap = 200

## Input files
minimap_alignment = 'SABE_1172_UNHESMSV_genotyping/SABE_1172_UNHESMSV_vs_GRCh38.sam.gz'
nrs_fasta_file = 'SABE_1172_UNHESMSV_merged/SABE1172_unmapped.fa'
contig_length_file = 'SABE_1172_UNHESMSV_ctg_lens_all.txt'

## Output files 
#If split has not been done, these will be created. If split was already run, the info will be pulled out of the files.
unmapped_output = "04_alignment_unmap.sam.gz"
partial_output = "04_alignment_partial.sam.gz"
#anchored_output = "04_alignment_anchored.sam.gz"
mapped_output = "04_alignment_map.sam.gz"

########################################################################################################################################

## Importing packages
import os
import gzip
import regex

# Check if file with contig lengths exist
if os.path.isfile(contig_length_file) == False:
    print("Creating file with contig lengths")
    os.system(f"awk '/^>/ {{if (seqlen){{print seqlen}}; printf $0\"\t\" ;seqlen=0;next; }} {{ seqlen += length($0)}}END{{print seqlen}}' {nrs_fasta_file} > {contig_length_file}")

## Separate alignment
# Check if outputfiles exist, if not, create them
dict_len = 0

if os.path.isfile(unmapped_output) == False or os.path.isfile(partial_output) == False or os.path.isfile(mapped_output) == False: # oros.path.isfile(anchored_output) == False:
    # Open alignment infile
    sam_in = gzip.open(minimap_alignment, "r")
    # Open alignment split outfiles
    unmapped = gzip.open(unmapped_output, "wb")
    mapped = gzip.open(mapped_output, "wb")
    partially = gzip.open(partial_output, "wb")
    #anchored = gzip.open(anchored_output, "wb")
    # Split alignments
    print("Splitting alignment")
    bmap = dict()
    for line in sam_in:
        #Skip header
        if line.startswith(b"@"):
            continue
        #Split input line
        ctg, flag, chrom, pos, mapq, cigar = line.split(b"\t")[:6]
        #Skip if secondary alignment
        if int(flag) & 2048:
            continue
        #Check divergence
        div, m = 1, None
        m = regex.search(b"de:f:(.+?)\s", line)
        if not m is None:
            div = float(m.group(1).decode())
        #Check cigar
        cig, m2 = None, None
        m2 = regex.findall("(\d{3,})[SI]", cigar.decode()) # Find all instances of SoftMask or Insertion in cigar
        if m2: 
            m2 = [float(i) for i in m2]
            cig = max(m2) #Work with the largest found instance
        # Split alignment
        if div > min_disc: #If divergence > min discordance => unmapped
            bmap[ctg.decode("utf-8")] = 1
            unmapped.write(line)
        elif (not cig is None) and (cig >= min_unmap): #If largest instance of S or I in cigar >= min unmapped sequence => partial
            bmap[ctg.decode("utf-8")] = 2
            partially.write(line)
            dict_len += 1
        else: #Else mapped
            bmap[ctg.decode("utf-8")] = 3
            mapped.write(line)
    # Close files
    sam_in.close()
    unmapped.close()
    mapped.close()
    partially.close()
    #
    print(" Done\n")
else: #If split files are in the folder, the partial mapping gets loaded in
    print("Alignment already split")
    print(" Creating alignment dictionary")
    bmap = dict()
    for line in gzip.open(partial_output, "rb"):
        ctg = line.split(b"\t")[0]
        bmap[ctg.decode("utf-8")] = 2
        dict_len += 1
    for line in gzip.open(unmapped_output, "rb"):
        ctg = line.split(b"\t")[0]
        bmap[ctg.decode("utf-8")] = 1
        dict_len += 1
    print(" Done\n")

## Binary map of alignment
print("Creating binary map of ctg alignment")
print(" Reading contig lengths")

bmap2 = dict()

for line in open(contig_length_file, "r"):
    ctg, length = line.split("\t")
    ctg = ctg.strip(">")
    if bmap.get(ctg) == 2: # if contig partially aligned
        bmap2[ctg] = "0" * int(length)
print("  Done")

##
print(" Processing partial alignment")

coordinates = dict()
sam_in = gzip.open(minimap_alignment, "r")

for line in sam_in:
    #Skip header
    if line.startswith(b"@"):
        continue
    #Split input line
    ctg, flag, chrom, pos, mapq, cigar = line.split(b"\t")[:6]
    ctg = ctg.decode("utf-8")
    #Skip if secondary alignment
    if int(flag) & 2048:
        continue
    #Check if contig in partial alignment contig dictionary; if not, skip
    if ctg not in bmap2:
        continue
    #Check divergence
    div = 1
    m3 = regex.search(b"de:f:(.+?)\s", line)
    if m3:
        div = float(m3.group(1).decode())
    if div > min_disc: # Skip unmapped contigs
        continue
    #Mark down coordinates of alignment
    coordinates[ctg] = f"{flag.decode('utf-8')}:{chrom.decode('utf-8')}:{pos.decode('utf-8')}:{mapq.decode('utf-8')}:{cigar.decode('utf-8')}"
    #Binary alignment map
    cur = 0
    alignlen = 0
    mapping = bmap2[ctg]
    #Reverse mapping if alignment reverse
    if int(flag) & 16:
        mapping = mapping[::-1]
    m4 = regex.findall(b"\d+\D", cigar) # Splits cigar into alignment instances
    for x in m4:
        alignlen = int(regex.search(b"\d+", x).group().decode()) # Extract alignment instance length
        if b"M" in x or b"I" in x: # Replace 0s in binary map with 1s
            mapping = mapping[:cur] + ("1" * alignlen) + mapping[cur + alignlen:]
        if b"D" not in x: # Deletion length skipped
            cur += alignlen
    #Reverse mapping if alignment reverse (put it back)
    if int(flag) & 16:
        mapping = mapping[::-1]
    #Save output to binary map
    bmap2[ctg] = mapping

sam_in.close()
print("  Done\n")

## Saving the output
print("Saving the output")

binaryoutput = gzip.open("04_binary_map.txt.gz", "wb")
i = 0
# j = 0

for contig, mapping in bmap2.items():
    pat1 = "^(0{" + str(min_unmap) + ",})"
    pat2 = "(0{" + str(min_unmap) + ",})$"
    if regex.search(pat1, mapping) or regex.search(pat2, mapping): #if there is long unmapped string at the beginning or end
        binaryoutput.write((f"{contig}\t{coordinates[contig]}\t{bmap2[contig]}\n").encode("utf-8"))
        i += 1
        # j += len(mapping)

##
print(f" Total partial: {i:,} contigs out of {dict_len:,}")
print("\nAll done, have a nice day.")
