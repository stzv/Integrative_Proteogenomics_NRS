from glob import glob
import gzip
import sys

files = ["02_mapped.sam.gz", "02_partial.sam.gz"]


list_output = open("03_Chm13v2_SABE_UNHESMSV_mapped_list.txt", "w")
list_output.write("#NRS,bp,mapping,flag,chromosome,position\n")

entries = []

for fil in files:
    print(f"Processing {fil}")
    ##
    alignment = fil[3:fil.find(".sam")]
    mpbs = 0
    ##
    f = gzip.open(fil,'rb')
    f = f.readlines()
    for line in f:
        nrs, flag, chr, pos = line.decode("utf-8").split("\t")[:4]
        nrs_length = (int(nrs.split("_")[-1]) - int(nrs.split("_")[-2]) + 1)/1000000
        entries.append(f'{nrs},{nrs_length},{alignment},{flag},{chr},{pos}')

if len(entries) == len(set(entries)):
    print("There are no duplicate IDs")
else:
    print("There are duplicate IDs, saving unique")

for item in set(entries):
    list_output.write(f"{item}\n")

list_output.close()
    
