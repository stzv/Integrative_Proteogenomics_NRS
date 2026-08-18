#!/usr/bin/env python

### ENVIRONMENT VARIABLES ####
min_cov = 7.5;
max_cov = 100;

## Input files
contig_length_file = 'SABE_1172_UNHESMSV_ctg_lens_all.txt'
unmapped_in = "04_alignment_unmap.sam.gz"
binary_map = "04_binary_map.txt.gz"
dark_cnt_loc = "brazil/20210114_genotyping/*_dark_cnt.txt.gz" 

## Output files 
freq_output = "05_frequency_results.txt"

########################################################################################################################################
## Importing packages
import gzip
import regex
from glob import glob
import time
from datetime import datetime

##
master_start = time.time()

log_fil = open(f"05_log_{datetime.now().strftime('%d%m%Y%H%M%S')}.txt", "w+")
print("Loading contigs' lengths")
start = time.time()

lens = dict()
for line in open(contig_length_file, "r"):
    ctg, length = line.split("\t")
    ctg = ctg.split(" ")[0]
    ctg = ctg.strip(">")
    lens[ctg] = int(length)

end = time.time()

print(f" Total: {len(lens):,} contigs\tElapsed time:{time.strftime('%Hh:%Mm:%Ss', time.gmtime(end-start))}")

log_fil.write("Contig lengths\n")
log_fil.write(f" Total: {len(lens):,} contigs\n")
log_fil.write(f"Elapsed time:{time.strftime('%Hh:%Mm:%Ss', time.gmtime(end-start))}\n\n")
## 
print("\nLoading unmapped contigs")
start = time.time()
total = 0
non38 = dict()

for line in gzip.open(unmapped_in, "r"):
    ctg = line.split(b"\t")[0]
    ctg = ctg.decode("utf-8")
    # Create dictionary with unmapped contigs with value 1
    non38[ctg] = 1
    total += lens[ctg]

end = time.time()

print(f" Total: {total/1000000:,.2f}Mbps in {len(non38):,} contigs\tElapsed time:{time.strftime('%Hh:%Mm:%Ss', time.gmtime(end-start))}")

log_fil.write("Unmapped contigs\n")
log_fil.write(f"Total: {total/1000000:,.2f}Mbps in {len(non38):,} contigs\n")
log_fil.write(f"Elapsed time:{time.strftime('%Hh:%Mm:%Ss', time.gmtime(end-start))}\n\n")

##
print("\nLoading partially mapped contigs")
start = time.time()

lens2 = dict()
counter = 0

for line in gzip.open(binary_map, "r"):
    counter += 1
    ctg, coordinates, bitstr = line.strip(b"\n").split(b"\t")
    ctg = ctg.decode("utf-8")
    if ctg == "k141_10000116":
        print(line)
    # Add to dictionary mapped contigs with value 2
    non38[ctg] = 2
    total += lens.get(ctg)
    # Checking for > 200 unmapped bps
    m1 = regex.match(b"^1+.*?(0{200,})$", bitstr) # Beginning of sequence
    m2 = regex.match(b"^(0{200,}).*?1+$", bitstr) # End of sequence
    if m1 is not None:
        alignlen = len(m1[1])
        lens2[ctg] = alignlen
    elif m2 is not None:
        alignlen = len(m2[1])
        lens2[ctg] = alignlen
    else:
        alignlen = lens.get(ctg)
        lens2[ctg] = alignlen
    if counter % 100000 == 0:
        print(f"{counter:,}")


end = time.time()

print(f" Total: {total/1000000:,.2f} Mbps in {len(non38):,} contigs\tElapsed time:{time.strftime('%Hh:%Mm:%Ss', time.gmtime(end-start))}")

log_fil.write("Partially unmapped contigs\n")
log_fil.write(f"Total: {total/1000000:,.2f} Mbps in {len(non38):,} contigs\n")
log_fil.write(f"Elapsed time:{time.strftime('%Hh:%Mm:%Ss', time.gmtime(end-start))}\n\n")
##
print("\nProcessing population frequency")
start = time.time()

table = dict() # dictionary of dictionaries - contains contig count for each sample
samples = []
total_len = 0

dark_cnt_files = glob(dark_cnt_loc) # list of all dark_cnt files
counter = 0

for fil in dark_cnt_files:
    counter += 1 #counter of how many files are processed
    new = 0 #counter of how many new contigs sample contains
    #Extract sample ID
    sample = regex.match(".*?(\d+)_dark_cnt.txt.gz", fil)[1]
    samples.append(sample)
    #Process contigs
    for line in gzip.open(fil, "r"):
        ctg, count = line.strip(b"\n").split(b"\t")
        ctg = ctg.decode("utf-8")
        count = int(count.decode("utf-8"))
        # Skip if contig is not in dictionary -> is mapped
        if ctg not in non38:
            continue
        #Create sample's ctg dictionary
        if ctg not in table:
            table[ctg] = {}
        if ('total' not in table[ctg]):
            table[ctg]['total'] = 0
        table[ctg][sample] = int(count) ## HERE OR ONLY IF CTG HAS PROPER COVERAGE??
        #Check for contig coverage
        if (count * 150 >= min_cov * lens[ctg]) and (count * 150 <= max_cov * lens[ctg]):
            new += 1
            table[ctg]['total'] += int(count)
            #Extract appropriate ctg length (unmapped or partial?)
            if ctg in lens2.keys():
                total_len += lens2[ctg]
            else:
                total_len += lens[ctg] #-> lens unmapped contigs
    print(f"  {sample} ({counter:,}/{len(dark_cnt_files):,})\tPartially unmapped contigs: {len(table):,}\tNew contigs: {new:,}\tTotal Mbps: {total_len/1000000:,.2f}")
    log_fil.write(f"{sample} ({counter:,}/{len(dark_cnt_files):,})\tPartially unmapped contigs: {len(table):,}\tNew contigs: {new:,}\tTotal Mbps: {total_len/1000000:,.2f}\n")

end = time.time()
log_fil.write(f"Elapsed time:{time.strftime('%Hh:%Mm:%Ss', time.gmtime(end-start))}\n\n")

## Finalizing the result
print("\nFinalizing the result")

#If contig is not in sample, it's count is 0
#Keep only contigs which have the value 'total' -> coverage within limits in at least one sample
table2 = table.copy() #Iterating through copy of dictionary -> avoid error 'dict has changed' when deleting contigs
counter = 0

print("\tRemoving contigs without proper coverage in at least one sample")
for contig, values in table2.items():
    if 'total' not in values or int(values['total']) == 0:
        del table[contig]

del table2

## Saving the result
print("\nCreating output")

frequency = open(freq_output, "w+")

total1 = 0
total2 = 0
bp1 = 0
bp2 = 0

#Print out numbers
for contig, vals in table.items():
    if non38[contig] ==1:
        bp1+= lens[contig]
        total1 += 1
    if non38[contig] == 2:
        bp2 += lens2[contig]
        total2 += 1
    frequency.write(f'{contig};{";".join(f"{k}:{v}" for k, v in vals.items())}\n')

master_end = time.time()

print(f"Fully unmapped: {total1:,} contigs, {bp1/1000000:,.2f} Mbps")
print(f"Partially unmapped: {total2:,} contigs, {bp2/1000000:,.2f} Mbps")

log_fil.write("\nResult\n")
log_fil.write(f"Fully unmapped: {total1:,} contigs, {bp1/1000000:,.2f} Mbps\n")
log_fil.write(f"Partially unmapped: {total2:,} contigs, {bp2/1000000:,.2f} Mbps\n")
log_fil.write(f"Total elapsed time:{time.strftime('%Hh:%Mm:%Ss', time.gmtime(master_end-master_start))}\n\n")
log_fil.close()

print("\nAll done, have a nice day")
