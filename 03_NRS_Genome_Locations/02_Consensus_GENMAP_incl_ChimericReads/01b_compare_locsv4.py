import sys
import collections
import os
from matplotlib_venn import venn3_unweighted, venn3_circles, venn3
from matplotlib import pyplot as plt

def Best_Anchor(anchors):
    all_anchors = anchors.split("|")[1:]
    anchor_count = len(all_anchors)
    all_anchors = sorted(all_anchors, key=lambda x: x.rsplit(':')[1])
    # Create & add new left and right anchor
    la = all_anchors[0]
    nla = ":".join([str(la.split(":")[0]), str(int(la.split(":")[1]) - 100), la.split(":")[2], la.split(":")[3]])
    all_anchors.insert(0, nla)
    ra = all_anchors[-1]
    nra = ":".join([str(ra.split(":")[0]), str(int(ra.split(":")[1]) - 100), ra.split(":")[2], ra.split(":")[3]])
    all_anchors.append(nra)
    # Best anchor location
    directions_read = [position.split(':')[2] for position in all_anchors]
    score_list = [(directions_read[:(i+1)].count('>') + directions_read[(i+1):].count('<')) for i in range(len(directions_read)-1)]
    best_pair = score_list.index(max(score_list))
    best_position = round((int(all_anchors[best_pair].split(':')[1]) + int(all_anchors[best_pair+1].split(':')[1])) / 2)
    # Best direction
    fwd_dir_count = directions_read.count(">")
    rvs_dir_count = directions_read.count("<")
    if fwd_dir_count >= rvs_dir_count:
        read_dir = "<"
    else:
        read_dir = ">"
    # Best strand
    direction_strand = [position.split(':')[3] for position in all_anchors]
    fwd_strand_count = direction_strand.count("+")
    rvs_strand_count = direction_strand.count("-")
    if fwd_strand_count >= rvs_strand_count:
        strand = "+"
    else:
        strand = "-"
    ##
    best_anchor = ":".join([str(la.split(":")[0]), str(best_position), str(anchor_count), read_dir, str(strand)])
    ##
    return best_anchor

# def Select_Best(anchors): # Based on read count
#     anchors_sorted = sorted(anchors, key=lambda x: x.rsplit(':')[2])
#     print(anchors_sorted)
    

########
Anchors_set, Anchors_mbp = set(), 0 ## Robustly mapped using Anchors
X10_set, X10_mbp = set(), 0 ## Robustly mapped using 10X
PMC_set, PMC_mbp = set(), 0 ## Mapped using PMC

NRS_dict = collections.defaultdict(dict)

chewed_loc_fil = "01b_SABE_1172_UNHESMSV_locations2.txt"

if os.path.isfile(chewed_loc_fil) == False:
    outfile = open(chewed_loc_fil, "w+")
    for line in open("01_SABE_1172_UNHESMSV_locations_table.tsv", "r"):
        if line.startswith("#") or line.startswith("NRS"):
            continue
        NRS, len_bp, PMC, X10, Anchors = line.strip("\n").split("\t")
        ##
        BestAnchor = "NA"
        if not Anchors == "NA":
            Anchors_set.add(NRS)
            Anchors_mbp += int(len_bp)
            BestAnchor = list()
            for anchor in Anchors.split(","):
                BestAnchor.append(Best_Anchor(anchor))
            BestAnchor = ",".join(BestAnchor)
        ##
        if not X10 == "NA":
            X10_set.add(NRS)
            X10_mbp += int(len_bp)
        ##
        if not PMC == "NA":
            PMC_set.add(NRS)
            PMC_mbp += int(len_bp)
        ##
        NRS_dict[NRS] = {"Length": len_bp, "PMC": PMC, "10X": X10.replace("chr", ""), "BestAnchor": BestAnchor}
        outfile.write("\t".join([NRS, len_bp, PMC, X10, BestAnchor]) + "\n")
else:
    for line in open(chewed_loc_fil, "r"):
        NRS, len_bp, PMC, X10, BestAnchor = line.strip("\n").split("\t")
        NRS_dict[NRS] = {"Length": len_bp, "PMC": PMC, "10X": X10.replace("chr", ""), "BestAnchor": BestAnchor}
        ##


unlocalized = set()
pmc_only = set()
anchor_only = set()
x10_only = set()

anchors_10x = set()
anchors_pmc = set()
x10_pmc = set()
anchor_10x_pmc = set()

### One location only
for nrs, info in NRS_dict.items():
    length, pmc, x10, anchors = info.values()
    #
    if pmc == "NA" and x10 == "NA" and "NA" in anchors: # No location
        unlocalized.add(nrs)
    elif not pmc == "NA" and x10 == "NA" and "NA" in anchors: # PMC only
        pmc_only.add(nrs)
    elif pmc == "NA" and x10 == "NA" and not "NA" in anchors: # Anchor only
        anchor_only.add(nrs)
    elif pmc == "NA" and not x10 == "NA" and "NA" in anchors: # 10X only
        x10_only.add(nrs)
    else:
        ### Anchor vs 10X
        if not x10 == "NA" and not "NA" in anchors:
            for anchor in anchors.split(","):
                x10_chr, x10_loc = x10.split(":")[:2]
                anchor_chr, anchor_loc, anchor_count = anchor.split(":")
                if anchor_chr == x10_chr and abs(int(x10_loc)*1000 - int(anchor_loc)) <= 100000:
                    anchors_10x.add(nrs)
                    continue # If multiple anchors, but at least 1 matche -> no need to check the rest
        ### Anchor vs PMC
        if not pmc == "NA" and not "NA" in anchors:
            for anchor in anchors.split(","):
                pmc_fl, pmc_chr, pmc_loc = pmc.split(":")[:3]
                anchor_chr, anchor_loc, anchor_count = anchor.split(":")
                if anchor_chr == pmc_chr and abs(int(pmc_loc) - int(anchor_loc)) <= 1000:
                    anchors_pmc.add(nrs)
                    continue # If multiple anchors, but at least 1 matche -> no need to check the rest
        ### PMC vs 10X
        if not pmc == "NA" and not x10 == "NA":
            pmc_fl, pmc_chr, pmc_loc = pmc.split(":")[:3]
            x10_chr, x10_loc = x10.split(":")[:2]
            if x10_chr == pmc_chr and abs(int(x10_loc)*1000 - int(pmc_loc)) <= 100000:
                x10_pmc.add(nrs)
        ### Anchor vs 10X vs PMC
        if not pmc == "NA" and not x10 == "NA" and not "NA" in anchors:
            for anchor in anchors.split(","):
                x10_chr, x10_loc = x10.split(":")[:2]
                anchor_chr, anchor_loc, anchor_count = anchor.split(":")
                pmc_fl, pmc_chr, pmc_loc = pmc.split(":")[:3]
                if x10_chr == pmc_chr and x10_chr == anchor_chr:
                    if abs(int(x10_loc)*1000 - int(anchor_loc)) <= 100000 and abs(int(pmc_loc) - int(anchor_loc)) <= 1000:
                        anchor_10x_pmc.add(nrs)
                        continue

all_bp = 0

for nrs in NRS_dict.keys():
    all_bp += int(NRS_dict[nrs]["Length"])
    check = "not good"
    for lst in pmc_only, anchor_only, x10_only, anchors_10x, anchors_pmc, x10_pmc, anchor_10x_pmc:
        if nrs in lst:
            check = "good"
    if check == "not good":
        unlocalized.add(nrs)

anchor_10x_nopmc = anchors_10x - anchor_10x_pmc 
anchor_pmc_no10x = anchors_pmc - anchor_10x_pmc 
x10_pmc_noanchor = x10_pmc - anchor_10x_pmc 

unlocalized_bp = round(sum([int(NRS_dict[nrs]["Length"]) for nrs in unlocalized])/ 1000000, 2)

pmc_only_bp = round(sum([int(NRS_dict[nrs]["Length"]) for nrs in pmc_only]) / 1000000, 2) # 001
anchor_only_bp = round(sum([int(NRS_dict[nrs]["Length"]) for nrs in anchor_only]) / 1000000, 2)#100
x10_only_bp = round(sum([int(NRS_dict[nrs]["Length"]) for nrs in x10_only]) / 1000000, 2) #011

anchor_10x_pmc_bp = round(sum([int(NRS_dict[nrs]["Length"]) for nrs in anchor_10x_pmc]) / 1000000, 2) # 111

anchor_10x_nopmc_bp = round(sum([int(NRS_dict[nrs]["Length"]) for nrs in anchor_10x_nopmc]) / 1000000, 2) #110
anchor_pmc_no10x_bp = round(sum([int(NRS_dict[nrs]["Length"]) for nrs in anchor_pmc_no10x]) / 1000000, 2) #101
x10_pmc_noanchor_bp = round(sum([int(NRS_dict[nrs]["Length"]) for nrs in x10_pmc_noanchor]) / 1000000, 2) #011


print("All NRS", all_bp)
print("Unlocalized", unlocalized_bp)
print(round(pmc_only_bp + anchor_only_bp + x10_only_bp + anchor_10x_pmc_bp + anchor_10x_nopmc_bp + anchor_pmc_no10x_bp + x10_pmc_noanchor_bp, 2))


# Percentages
all_bp = all_bp/1000000
anchor_only_perc = round((anchor_only_bp * 100.0) / all_bp, 2)
x10_only_perc = round((x10_only_bp * 100.0) / all_bp, 2)
anchors_10x_mbp_perc = round((anchor_10x_nopmc_bp * 100.0) / all_bp, 2)
pmc_only_perc = round((pmc_only_bp * 100.0) / all_bp, 2)
anchors_pmc_mbps_perc = round((anchor_pmc_no10x_bp * 100.0) / all_bp, 2)
x10_pmc_mbp_perc = round((x10_pmc_noanchor_bp * 100.0) / all_bp, 2)
anchor_10x_pmc_mpbs_perc = round((anchor_10x_pmc_bp * 100.0) / all_bp, 2)


fig, ax = plt.subplots()

v = venn3_unweighted(subsets=(anchor_only_bp, x10_only_bp, anchor_10x_nopmc_bp, pmc_only_bp, anchor_pmc_no10x_bp, x10_pmc_noanchor_bp, anchor_10x_pmc_bp), 
                    set_labels=('Anchors', '10X', 'PMC'), 
                    set_colors=("orange", "blue", "red"), alpha=0.7)

v.get_label_by_id("100").set_text(str(anchor_only_bp) + "Mbp\n" + str(anchor_only_perc) + "%")
v.get_label_by_id("110").set_text(str(anchor_10x_nopmc_bp) + "Mbp\n" + str(anchors_10x_mbp_perc) + "%")
v.get_label_by_id("010").set_text(str(x10_only_bp) + "Mbp\n" + str(x10_only_perc) + "%")
v.get_label_by_id("101").set_text(str(anchor_pmc_no10x_bp) + "Mbp\n" + str(anchors_pmc_mbps_perc) + "%")
v.get_label_by_id("111").set_text(str(anchor_10x_pmc_bp) + "Mbp\n" + str(anchor_10x_pmc_mpbs_perc) + "%")
v.get_label_by_id("011").set_text(str(x10_pmc_noanchor_bp) + "Mbp\n" + str(x10_pmc_mbp_perc) + "%")
v.get_label_by_id("001").set_text(str(pmc_only_bp) + "Mbp\n" + str(pmc_only_perc) + "%")

plt.suptitle("SABE_1172_UNHESMSV")
plt.savefig("01b_SABE_1172_UNHESMSV_Venn_Locations_NRSsumMbps2.png")