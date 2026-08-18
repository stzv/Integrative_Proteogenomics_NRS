from glob import glob
import sys

#infiles = glob("/mnt/fedot21/brazil/20210622_STAR/results/*/*Chimeric.out.junction")

infiles = list()
for line in open("01_SABE_UNHESMSV_1172_RNASeq_infile_list.txt", "r"):
    infil = line.strip("\n").replace("Log.final.out", "Chimeric.out.junction")
    infiles.append(infil)

outfile = open("05_RNASeq_chimeric_junctions_merged_binned.txt", "w+")
counter = 0

print(" Total RNA Files:", len(infiles))

for infile in infiles:
    #rna_id = infile[:infile.find("Chimeric.out.junction")]
    counter += 1
    if counter % 1000 == 0:
        print(f"{counter}/{len(infiles)}")
    rna_sample = infile.split("/")[-1].strip("Chimeric.out.junction")
    cj_dict = dict()
    for line in open(infile, "r"):
        d_chrom, d_position = line.split("\t")[:2]
        donor = d_chrom + ":" + d_position
        a_chrom, a_position = line.split("\t")[3:5]
        acceptor = a_chrom + ":" + a_position
        read_name = line.split("\t")[9]
        ##
        if "k141_" in d_chrom or "k141_" in a_chrom:
            if ("k141" in acceptor) and ("k141" in donor): # Skip NRS to NRS Chimera
                continue
            elif "k141" in donor:
                nrs = donor.split(":")[0]
                nrs_loc = donor
                hg38_loc = acceptor
            elif "k141" in acceptor:
                nrs = acceptor.split(":")[0]
                nrs_loc = acceptor
                hg38_loc = donor
            ## Bin the CJ
            # Prepare values
            chrom = hg38_loc.split(":")[0]
            pos = hg38_loc.split(":")[1]
            if nrs not in cj_dict.keys():
                cj_dict[nrs] = {}
            # Create a bucket
            if hg38_loc.startswith("HLA"):
                bucket = hg38_loc
                cj_dict[nrs].setdefault(bucket, 0)
            else:
                if len(pos) < 4:
                    bucket = 500
                    cj_dict[nrs].setdefault(f"{chrom}:{bucket}-{bucket+500}", 0) # bucket with higher range
                    cj_dict[nrs].setdefault(f"{chrom}:{bucket-500}-{bucket}", 0) # bucket with lower range
                else:
                    bucket = round(int(pos),-3) # Create the bucket multiple of 1000 
                    cj_dict[nrs].setdefault(f"{chrom}:{bucket}-{bucket+500}", 0) # bucket with higher range
                    cj_dict[nrs].setdefault(f"{chrom}:{bucket-500}-{bucket}", 0) # bucket with lower range
            ## Check into which bucket CJ belongs
            for b in cj_dict[nrs]:
                if hg38_loc.startswith("HLA"):
                    if hg38_loc == b: #HLA hits
                        cj_dict[nrs][b] += 1
                else:
                    if (chrom == b.split(":")[0]) and (int(b.split(":")[1].split("-")[0]) <= int(pos)) and (int(b.split(":")[1].split("-")[1]) >= int(pos)):
                        cj_dict[nrs][b] += 1
    ####
    for nrs in cj_dict:
        bin_list = [f"{b}|{cj_dict[nrs][b]}" for b in cj_dict[nrs]]
        outfile.write(f"{nrs}\t{','.join(bin_list)}\t{rna_sample}\n")

import subprocess
subprocess.call("email_stepanka.pl merge chimeric junctions", shell = True)

print("\nAll done, have a nice day!")