import gzip
import collections
import sys

alignmentfil = '00_UN_vs_UNHESMSV.m6'

####
print("Load in UN library")
NRS_UN_dict = collections.defaultdict(int)
for line in open("Sao_Paulo_dark_freeze2.fa", "r"):
    if line.startswith(">"):
        nrs = line.strip("\n").replace(">", "").split()[0]
        length = int(nrs.split("_")[-1]) - int(nrs.split("_")[-2]) + 1
        NRS_UN_dict[nrs] = length


####
print("Load in alignment")
UN_BLAST_dict = collections.defaultdict(list)
all_nrs_in_blast, kept_in_blast = 0, 0
count = 0

for line in open(alignmentfil, "r"):
    qseqid, sseqid, pident, length, qstart, qend, sstart, send, evalue = line.split("\t")
    qseqid_len = int(qseqid.split("_")[-1]) - int(qseqid.split("_")[-2]) + 1
    if (int(length)/qseqid_len) < 0.8:
        continue
    UN_BLAST_dict[qseqid].append([sseqid, float(pident), int(length), int(qstart), int(qend), sstart, send, qseqid_len])
    ##
    count += 1

outfil1 = open("01a_UN_notAlignedByBlast.txt", "w+")
UN_notBLAST_count, UN_notBLAST_mbp = set(), 0
for UN_nrs, length in NRS_UN_dict.items():
    if UN_nrs not in UN_BLAST_dict.keys():
        UN_notBLAST_mbp += length
        outfil1.write(f"{UN_nrs}\t{length}\n")


####
print("Choose best alignment")
NRS_UNHESMSV_dict = collections.defaultdict(list)

for UN_nrs, info in UN_BLAST_dict.items():
    best_UN_IM, best_UN_alignment = 0, ""
    nrs_len = info[0][0]
    ##
    best_UN_IM = max([float(sublist[2]) for sublist in info])
    for lst in info:
        if best_UN_IM == float(lst[2]):
            best_UN_alignment = lst
    # ##
    unhesmsv_nrs = best_UN_alignment[0]
    unhesmsv_info = [UN_nrs] + best_UN_alignment[1:]
    NRS_UNHESMSV_dict[unhesmsv_nrs].append(unhesmsv_info)
    ##

####
print("Count alignment types")
count_expanded, count_expanded_mbp = set(), 0
count_carried, count_carried_mbp = set(), 0
un_merged_count, un_merged_mbp, un_lost_mbp = set(), 0, 0

for nrs, info in NRS_UNHESMSV_dict.items():
    # Get NRS length
    nrs_len = int(nrs.split("_")[-1]) - int(nrs.split("_")[-2]) + 1
    #print(nrs)
    # Expands only one UN contig
    if len(info) == 1: 
        #print(info)
        count_expanded.add(nrs)
        nrs_un, im, al_len, un_start, un_end, unhesmsv_start, unhesmsv_end, nrs_un_len = info[0]
        alignment = abs(un_start - un_end + 1)
        # Add how much bp the UNHESMSV contig adds extra to the UN contig
        count_expanded_mbp += nrs_len - alignment
        count_carried_mbp += alignment
        un_merged_mbp += alignment
        un_lost_mbp += nrs_un_len - abs(un_start - un_end + 1) # part of UN which does not align to the UNHESMSV
        #print(un_lost_mbp)
    # Merges and expands multiple UN contigs
    elif len(info) > 1:
        count_carried.add(nrs)
        cov_list = list()
        # For each UN contig, keep which portion of UNHESMSV contig is covered
        #print(nrs)
        for i in info:
            #print(i)
            nrs_un, im, al_len, un_start, un_end, unhesmsv_start, unhesmsv_end, nrs_un_len = i
            alignment = abs(un_start - un_end + 1)
            #print(i)
            un_merged_mbp += alignment
            un_lost_mbp += nrs_un_len - abs(un_start - un_end + 1) # part of UN which does not align to the UNHESMSV
            #un_merged_count.add(nrs_un)
            cov_list = cov_list + [int(unhesmsv_start), int(unhesmsv_end)]
            #print(un_lost_mbp)
        # Get the start and end of the UNHESMSV coverage
        cov_start = min(cov_list)
        cov_end = max(cov_list)
        # Calculate how much of UNHESMSV contig is covered
        nrs_cov = cov_end - cov_start + 1
        # Calculate how much is the UNHESMSV contig expanding
        expanded = nrs_len - nrs_cov
        count_expanded_mbp += nrs_len - nrs_cov
        count_carried_mbp += nrs_cov
    else:
        print("Something weird is happening with contig", nrs)
        sys.exit()

####
print("Load in UNHESMSV library")
count_total, total_mbp = set(), 0
count_new, count_new_mbp = set(), 0

for line in open("../01_NRS_Assembly/SABE_1172_UNHESMSV_NRS_dark_freeze_final.fa", "r"):
    if line.startswith(">"):
        nrs = line.split()[0].replace(">", "")
        total_mbp += int(nrs.split("_")[-1]) - int(nrs.split("_")[-2]) + 1
        count_total.add(nrs)
        if nrs not in count_expanded and nrs not in count_carried:
            count_new.add(nrs)
            count_new_mbp += int(nrs.split("_")[-1]) - int(nrs.split("_")[-2]) + 1

#print(count_expanded)
print("\n")

print("\tTotal", len(count_total), "NRS", round(total_mbp/1000000, 2), "Mbp")
print("\tExpanded", len(count_expanded), "NRS", round(count_expanded_mbp/1000000, 2), "Mbp")
print("\tCarried over", len(count_carried), "NRS", round(count_carried_mbp/1000000, 2), "Mbp")
print("\tNew", len(count_new), "NRS", round(count_new_mbp/1000000, 2), "Mbp")

print("Mbp not accounted for", (round((total_mbp - (count_expanded_mbp + count_carried_mbp + count_new_mbp))/1000000, 2)))

print("\n")

un_total = sum(list(NRS_UN_dict.values()))
print("\tTotal UN library", len(NRS_UN_dict.keys()), "NRS", round(un_total/1000000, 2), "Mbp")
print("\tUN NRS not aligned by BLAST (<0.8)", round(UN_notBLAST_mbp/1000000, 2), "Mbp")
print("\tUN merged together/carried over", round(un_merged_mbp/1000000, 2), "Mbp")
print("\tUN lost", round(un_lost_mbp/1000000, 2), "Mbp")

print("Mbp not accounted for", (round((un_total - (UN_notBLAST_mbp + un_merged_mbp + un_lost_mbp))/1000000, 2)))

# accounted_for = len(count_expanded) + len(count_carried) + len(count_new)
# accounted_for_mbp = count_expanded_mbp + count_carried_mbp + count_new_mbp

# print("\n")
# print("Sum",accounted_for, "NRS", round(accounted_for_mbp/1000000, 2), "Mbp")
# print("Leftover", len(count_total) - accounted_for, "NRS", round((total_mbp - accounted_for_mbp)/1000000, 2), "Mbp")
# #print("\nMerged", len(count_carried))

