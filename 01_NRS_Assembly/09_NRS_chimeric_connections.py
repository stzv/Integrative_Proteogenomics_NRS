import subprocess
import os
import collections
import regex

### ENVIRONMENT VARIABLES ####
# min_disc = 0.05
# min_unmap = 200

## Input files
minimap_alignment = 'SABE_1172_UNHESMSV_vs_GRCh38.sam.gz'
freeze_fasta = "SABE_1172_UNHESMSV_NRS_dark_freeze_final.fa"
NRS_ID_file = "NRS_ID_list.txt"

align_PU = "04_alignment_partial.sam.gz" #04_alignment_unmapped.sam.gz
NRS_ID_file_PU_all = "09_NRS_ID_list_partial_unmapped_all.txt"
NRS_ID_file_PU = "09_NRS_ID_list_partial.txt"
Alignments_file = "SABE_1172_UNHESMSV_vs_GRCh38.sam.gz"

SuppAligments_file = "09_SuppSec_Alignments.txt"
PrimSuppAlignment_file = "09_PrimSupp_Alignments.txt"

####
## In 04_split alignment, the ctgs were split based on whether they are mapped, unmapped, or partially mapped
## Unmapped: divergence >= 0.05 (< 95% IM)
## Partially: divergence < 0.5, largest instance of M or I >= 200bp
## Mapped: divergence < 0.5, I or M portion < 200bp
##
## The minimap alignments were filtered out for flag 2048 = supplementary alignment (chimeric)
## Some partial alignments could have secondary alignment
##
## Example NRS k141_19408630
## k141_19408630   0       2       218581274       60      626M614S
## k141_19408630   2048    8       25196971        60      580H 68M 1D 71M 21D 433M 4I 11M 7D 73M
##
## In genome mapping, the Chr8 portion of the NRS maps with 94% IM, forward strand, exactly according to the CIGAR (H = hard clipping)
## Sum of matches in the cigar = 656M, the length of the "NRS" portion = 614, the whole NRS length = 1240bp
## 1240 - 580H = 660 bp ; 656/660 = 99% completeness (not computing the penalty for I or D)
## Could be considered chimeric alignment => could be GRIP, assembly artifact, or structural variant
## As such, we cannot conlude what exactly it is, therefore we have to put them aside from the final results

## !!!!!
## Filtering should be applied on NRS, whose secondary alignment aligns the Soft/Hard-masked portion of the primary alignment

## 
print("Get the NRS which were marked as partial mapped")
if os.path.isfile(NRS_ID_file_PU) == False:
    command = f"zcat {align_PU} | samtools view -t SABE_1172_UNHESMSV_alignment_HEADER.txt -F 2048 -q 60 | grep -f {NRS_ID_file} | grep 'SA:Z:' | awk '{{print $1,$2,$3,$4,$5,$6,$21}}' > {NRS_ID_file_PU}"
    print(command)
    subprocess.call(command, shell = True)

##
print("Load the contig lengths")
NRS_len_dict = collections.defaultdict(dict)
for line in open("SABE_1172_UNHESMSV_ctg_lens_all.txt", "r"):
    nrs, lenght = line.strip("\n").split("\t")
    nrs = "_".join(nrs.replace(">", "").split(" ")[0].split("_")[:2])
    NRS_len_dict[nrs] = int(lenght)

##
print("Process the alignments")

def BinaryMap(cigar, cur, flag, mapping):
    m4 = regex.findall("\d+\D", cigar) # Splits cigar into alignment instances
    for x in m4:
        alignlen = int(regex.search("\d+", x).group()) # Extract alignment instance length
        if "M" in x or "I" in x: # Replace 0s in binary map with 1s
            mapping = mapping[:cur] + ("1" * alignlen) + mapping[cur + alignlen:]
        if "D" not in x: # Deletion length skipped
            cur += alignlen
    ##
    return mapping

def ProcessNRS(line):
    nrs, flag, chrom, pos, mapq, cigar, secondary_alignment = line.strip("\n").split()
    sec_als = secondary_alignment.split(";")
    length = NRS_len_dict[nrs]
    ## Prepare variables
    bmap = "0" * int(length)
    cur = 0
    ## First create binary map of the primary alignment
    bmap = BinaryMap(cigar, cur, flag, bmap)
    ## Add the mappings from the secondary alignments
    for sec_al in sec_als:
        if not sec_al: continue
        cur = 0
        s_chr, s_pos, s_flag, s_cigar, s_mapq, s_num = sec_al.replace("SA:Z:", "").split(",")
        bmap = BinaryMap(s_cigar, cur, s_flag, bmap)
    ## Check if there is still >= 200bp of unmapped sequence
    m1 = regex.match("^1+.*?(0{200,})$", bmap) # Beginning of sequence
    m2 = regex.match("^(0{200,}).*?1+$", bmap) # End of sequence
    if not m1 and not m2: # no more unmapped sequence
        outline = "\t".join([nrs, chrom, pos, cigar, secondary_alignment])
        outfile.write(outline + "\n")

NRS_Al_dict = collections.defaultdict(list)
#counter = 0
progress_counter = 0
outfile = open("09_NRS_SecMapping_resolved.txt", "w")
for line in open(NRS_ID_file_PU, "r"):
    if not 'SA:Z:' in line: continue # If it happens that there is a line without secondary alignments
    try:
        ProcessNRS(line)
    except Exception:
        print("ERROR")
        print(line)
        continue
    ##
    progress_counter += 1
    if progress_counter % 1000 == 0:
        print(progress_counter)

