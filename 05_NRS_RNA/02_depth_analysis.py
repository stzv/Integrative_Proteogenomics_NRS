from glob import glob
import subprocess
import os
import more_itertools as mit
import re
import pysam


## ENVIROMENT VARIABLES
## Minimum continuous coverage
minconcov = 10
## Minimum aligned reads
minreads = 10
## Min quality for depth calculation
Q = 20

## INPUT FILES


print("Sort & Filter BAM files")
files = "01_SABE_UNHESMSV_1172_RNASeq_infile_list.txt"
bamfiles = "02_SABE_1172_UNHESMSV_RNASeq_bamfiles.txt"

rna_samples = list()
counter = 0 

with open(bamfiles, "w+") as bamout:
    for line in open(files, "r"):
        counter += 1
        ## Extract RNAID
        rid = line.strip("\n").split("/")[-1].strip("Log.final.out")
        rna_samples.append(rid)
        ### Prepare file names
        bamfil = line.strip("\n").replace("Log.final.out", ".bam")
        sortedbamfil = line.strip("\n").replace("Log.final.out", "_sorted.bam").split("/")[-1]
        filteredbamfil = sortedbamfil.replace("sorted", "sorted_filtered")
        bamout.write(f'{filteredbamfil}\n')
        ### Sort bam file
        if not os.path.exists(sortedbamfil) and not os.path.exists(filteredbamfil):
            ## Sort bam file
            print(counter, sortedbamfil)
            subprocess.call(f"samtools sort {bamfil} > {sortedbamfil}", shell = True)
        ### Remove reads with >= 10SM ends
        if not os.path.exists(filteredbamfil):
            print(counter, filteredbamfil)
            samfile = pysam.AlignmentFile(sortedbamfil, "rb") # Infile
            header = samfile.header
            outfile = pysam.AlignmentFile(filteredbamfil, "wb", header = header)# Outfile
            ## Read the alignment
            for read in samfile:
                ## Start with clean slate
                SM_match = ""
                save = "YES" # If cigar starts or ends with Soft Masked, will change to no and will be ommited
                ## Check if Soft Masked in CIGAR
                if "S" in read.cigarstring:
                    SM_match = re.findall("(\d{2,})[S]", read.cigarstring) # 2 digits => >=10S
                    # Check if cigar starts or ends with softmasked
                    for match in SM_match:
                        m = match + "S"
                        if read.cigarstring.startswith(m) or read.cigarstring.endswith(m):
                            save = "NO"
                ## Output read into new bam file
                if save == "YES":
                    outfile.write(read)
        ## Remove sorted file
        if os.path.exists(sortedbamfil):
            print(f"Remove tmp file {sortedbamfil}")
            subprocess.call(f"rm {sortedbamfil}", shell = True)


print("Extract reads count")
coverage, continuity, contin, maxcoverage = dict(), dict(), dict(), dict()
counter_rd = 0

for line in open(bamfiles, "r"):
    # Run samtools depth
    infile = line.strip("\n")
    depthout = infile.replace("_sorted_filtered.bam", "_depth.txt")
    subprocess.call(f"samtools depth {infile} -Q {Q} > {depthout}", shell = True)
    # Extract RNA Sample ID & ascenssion number    
    rnaid = infile.replace("_sorted_filtered.bam", "")
    # For table creation
    ## Process depth file
    # Extract max # of reads aligned to 1bp in NRS -> coverage dictionary
    # Extract how many continuous bp is covered in NRS -> continuity
    for entry in open(depthout, "r"):
        ## Extract ctg ID & add to dictionaries
        ctg = entry.split("\t")[0]
        if ctg not in coverage.keys():
            coverage[ctg] = {}
        if ctg not in continuity.keys():
            continuity[ctg] = {}
        if ctg not in contin.keys():
            contin[ctg] = {}
        if ctg not in maxcoverage.keys():
            maxcoverage[ctg] = {}
        ## Add RNAID to dictionaries
        if rnaid not in coverage[ctg]:
            coverage[ctg][rnaid] = []
        if rnaid not in continuity[ctg]:
            continuity[ctg][rnaid] = []
        if rnaid not in contin[ctg]:
            contin[ctg][rnaid] = "FALSE" # Consider continuity false by default
        #if rnaid not in maxcoverage[ctg]:
        #    maxcoverage[ctg][rnaid] = 0
        ## Save the numbers
        readcount = int(entry.split("\t")[2])
        pos = int(entry.split("\t")[1])
        ## Evaluate coverage and continuity
        if readcount >= minreads: # No of aligned reads to basepair > minreads
            coverage[ctg][rnaid].append(readcount) # Append the depth to dictionary
            continuity[ctg][rnaid].append(pos) # Append the covered position
        ## Evaluate whether ctg's continuity in RNAID is >= min continuity count
        cont_list = [list(group) for group in mit.consecutive_groups(continuity[ctg][rnaid])] # Split positions into consecutive groups
        for group in cont_list:
            if len(group) >= minconcov:
                contin[ctg][rnaid] = "TRUE"
                ## Save the highest read depth from continuously covered sequence
                contin_index = continuity[ctg][rnaid].index(group[-1]) # Get index of new position in sequence
                maxreads = coverage[ctg][rnaid][contin_index] # Fetch the readcount of this position
                if maxreads > maxcoverage[ctg].get(rnaid, 0): #maxcoverage[ctg][rnaid]: # If readcount higher than previous, save as new max
                    maxcoverage[ctg][rnaid] = maxreads
    ## Counter
    #counter_rd += 1
    #if counter_rd == 10:
    #    break       
    ## Remove temporary file with depths
    subprocess.call(f"rm {depthout}", shell = True)


#print(maxcoverage)
# Free up memory
del contin, coverage, continuity

print("Saving to a file")
outfile = open("02_SABE_1172_UNHESMSV_NRS_RNA_reads_count.txt", "w+")
delim = "\t"
outfile.write(f"#RNA_list\t{delim.join(rna_samples)}\n")


for ctg, nested in maxcoverage.items():
    ## For contigs not covered by RNASample, add 0 as maxcoverage
    #for rnaid in rna_samples:
    #    if rnaid not in values:
    #        maxcoverage[ctg][rnaid] = 0
    ## Sort coverage according to RNA list in order to print out table
    #maxcoverage[ctg] = {k: maxcoverage[ctg][k] for k in rna_samples}
    ## Print out
    #counts = delim.join(str(v) for c, v in maxcoverage[ctg].items())
    #outfile.wite(f"{ctg}\t{counts}\n")
    outfile.write(f"{ctg}")
    for subkey, value in nested.items():
        outfile.write(f"\t{subkey}:{value}")
    outfile.write("\n")

subprocess.call("email_stepanka.pl depth analysis", shell = True)
print("All done, have a nice day!")