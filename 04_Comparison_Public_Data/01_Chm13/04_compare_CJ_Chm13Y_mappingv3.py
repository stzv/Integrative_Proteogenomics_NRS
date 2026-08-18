import sys
import collections
import os
#from matplotlib_venn import venn3_unweighted, venn3_circles, venn3
from matplotlib import pyplot as plt
import pandas as pd
from upsetplot import from_memberships
from upsetplot import UpSet
import upsetplot
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import sys

#####
print("Compare locations #1")
RNASeq_dict = collections.defaultdict(dict)

for line in open("03_NRS_Genome_Locations/02_Consensus_GENMAP_incl_ChimericReads/06_SABE_1172_UNHESMSV_NRS_best_CJ.txt", "r"):
    #k141_49970037_1_436	MT:11000-11500|82305
    nrs, locs = line.strip("\n").split("\t")
    nrs_short = "_".join(nrs.split("_")[:2])
    if locs.startswith("HLA"):
        continue
    RNASeq_dict[nrs_short] = locs

####
CHM13_dict = collections.defaultdict(list)
for line in open("03_Chm13v2_SABE_UNHESMSV_mapped_list.txt", "r"):
    if line.startswith("#"):
        continue
    nrs, nrs_len, mapping, flag, chr, pos = line.strip("\n").split(",")
    nrs_short = "_".join(nrs.split("_")[:2])
    CHM13_dict[nrs_short].append(":".join([chr.replace("chr", ""), pos]))


###### Upsetplot - only whether they have location or not

NRS_dict = collections.defaultdict(int)
LOC_dict = collections.defaultdict(set)

for line in open("03_NRS_Genome_Locations/02_Consensus_GENMAP_incl_ChimericReads/01b_SABE_1172_UNHESMSV_locations2.txt", "r"):
    nrs, bp, pmc, x10, anchors = line.strip("\n").replace("chr", "").split("\t")
    NRS_dict[nrs] = int(bp)/1000000
    # Get the locations from other sources
    rna = "NA"
    if nrs in RNASeq_dict.keys():
        rna = RNASeq_dict[nrs]
    chm13 = "NA"
    if nrs in CHM13_dict.keys():
        chm13 = CHM13_dict[nrs]
    # Add the locations to a dictionary
    if not pmc == "NA":
        LOC_dict["PMC"].add(nrs)
    if not anchors == "NA":
        LOC_dict["Anchors"].add(nrs)
    if not rna == "NA":
        LOC_dict["RNASeq"].add(nrs)
    if not x10 == "NA":
        LOC_dict["10X"].add(nrs)
    if not chm13 == "NA":
        LOC_dict["Chm13"].add(nrs)
    ## No location at all
    if all(x  == "NA" for x in (pmc, x10, anchors, rna, chm13)):
        LOC_dict["Unloc"].add(nrs)


### preparation for upsetplot

nrs_df = pd.DataFrame(list(NRS_dict.items()))
nrs_df.columns = ["nrs", "length"]

chm13_df = pd.DataFrame(list(LOC_dict["Chm13"]))
chm13_df.columns = ["nrs"]
chm13_df["length"] = chm13_df['nrs'].map(NRS_dict)
print("Total NRS mapped to CHM13:", len(chm13_df))
print("Fraction of total NRS mapped to CHM13:",
      round(len(chm13_df) / len(NRS_dict) * 100, 2), "%")


anch_df = pd.DataFrame(list(LOC_dict["Anchors"]))
anch_df.columns = ["nrs"]
anch_df["length"] = anch_df['nrs'].map(NRS_dict)
anch_df = anch_df[(anch_df["nrs"].isin(chm13_df["nrs"]))]

pmc_df = pd.DataFrame(list(LOC_dict["PMC"]))
pmc_df.columns = ["nrs"]
pmc_df["length"] = pmc_df['nrs'].map(NRS_dict)
pmc_df = pmc_df[(pmc_df["nrs"].isin(chm13_df["nrs"]))]

x10_df = pd.DataFrame(list(LOC_dict["10X"]))
x10_df.columns = ["nrs"]
x10_df["length"] = x10_df['nrs'].map(NRS_dict)
x10_df = x10_df[(x10_df["nrs"].isin(chm13_df["nrs"]))]

rna_df = pd.DataFrame(list(LOC_dict["RNASeq"]))
rna_df.columns = ["nrs"]
rna_df["length"] = rna_df['nrs'].map(NRS_dict)
rna_df = rna_df[(rna_df["nrs"].isin(chm13_df["nrs"]))]

noloc_df = pd.DataFrame(list(LOC_dict["Unloc"]))
noloc_df.columns = ["nrs"]
noloc_df["length"] = rna_df['nrs'].map(NRS_dict)


###
print("Get NRS mapped by Chm13 but no other GENMAP")
chm13_only_df = chm13_df[~(chm13_df["nrs"].isin(anch_df["nrs"]))
                       & ~(chm13_df["nrs"].isin(pmc_df["nrs"]))
                       & ~(chm13_df["nrs"].isin(x10_df["nrs"]))
                       & ~(chm13_df["nrs"].isin(rna_df["nrs"]))
                       & ~( chm13_df["nrs"].isin(noloc_df["nrs"]))
                        ]

nrs_chm13_only = set(chm13_only_df["nrs"].to_list())
#chm13_only_df.to_csv("04_NRS_only_Chm13.txt")
outfile = open("04_NRS_only_Chm13.txt", "w+")

for nrs in NRS_dict.keys():
    if nrs in nrs_chm13_only:
        outfile.write("\t".join([nrs, ",".join(CHM13_dict[nrs])]) + "\n")

# Draw the uspetplot
print("Create upsetplot #1")

content2 = {"Total NRS": [nrs_df["nrs"]],
        #"Unloc NRS": [noloc_df["nrs"]],
        "Anchors": [anch_df["nrs"]],
        "PMC": [pmc_df["nrs"]],
        "10X": [x10_df["nrs"]],
        "RNASeq": [rna_df["nrs"]],
        "Chm13":   [chm13_df["nrs"]],
        }

data2 = upsetplot.from_contents(content2)
#data2 = upsetplot.query(data2, present = ["Chm13"]).data


#
print("Draw upsetplot #1")
fig = plt.figure(figsize = (6, 2)) # width, height
ax = fig.add_subplot()

Upset_Intersec = UpSet(data2,
                        sort_by = "cardinality", 
                        subset_size = "count", 
                        show_counts = True,
                        #show_percentages = True,
                        #min_subset_size = 6000, # > 1% of NRS
                        element_size = 30
                        )


Upset_Intersec.style_subsets(
        present = ["Chm13"],#, "Total NRS"],
        absent = ["Anchors", "PMC", "10X", "RNASeq"],
        facecolor = "chocolate"
        )

Upset_Intersec.style_subsets(
        present = ["Total NRS"],#, "Total NRS"],
        absent = ["Chm13", "Anchors", "PMC", "10X", "RNASeq"],
        facecolor = "green"
        )

params = {'font.size': 9}
with plt.rc_context(params):
    Upset_Intersec.plot()

plt.text(0, 1.5, "Mapped NRS #", 
        fontsize = 10, 
        transform = ax.transAxes,
        backgroundcolor = 'white'
        )

#plt.suptitle("SABE_1172_UNHESMSV")
plt.savefig("04_SABE_UNHESMSV_Genmap_RNASeq_Chm13Y_Upset_LocationsCount.png")


# ###### Upsetplot - match the locations
# print("Compare locations #2")

# NRS_dict = collections.defaultdict(int)
# LOC_dict = collections.defaultdict(set)

# MATCH_dict = list()

# for line in open("SABE_1172_UNHESMSV_genomemapping/01b_SABE_1172_UNHESMSV_locations2.txt", "r"):
#     ## Prepare variables
#     nrs, bp, pmc, x10, anchors = line.strip("\n").replace("chr", "").split("\t")
#     NRS_dict[nrs] = int(bp)/1000000
#     ##
#     rna = "NA"
#     if nrs in RNASeq_dict.keys():
#         rna = RNASeq_dict[nrs]
#     chm13 = "NA"
#     if nrs in CHM13_dict.keys():
#         chm13 = CHM13_dict[nrs]
#     ##
#     ## Match the locations
#     # Anchor
#     if not "NA" in anchors:
#         reported_match = set()
#         for anchor in anchors.split(","):
#             anchor_chr, anchor_loc = anchor.split(":")[:2]
#             # Vs 10X
#             if not x10 == "NA":
#                 x10_chr, x10_loc = x10.split(":")[:2]            
#                 if anchor_chr == x10_chr and abs(int(anchor_loc) - (int(x10_loc)*1000)) <= 100000:
#                     reported_match.update(["Anchor", "10X"])
#             # Vs PMC
#             if not pmc == "NA":
#                 pmc_fl, pmc_chr, pmc_loc = pmc.split(":")[:3]
#                 if anchor_chr == pmc_chr and abs(int(anchor_loc) - int(pmc_loc)) <= 1000:
#                     reported_match.update(["Anchor", "PMC"])
#             # Vs RNASeq
#             if not rna == "NA":
#                 for r in rna.split(","):
#                     rna_chr = r.split("|")[0].split(":")[0]
#                     rna_start, rna_end = r.split("|")[0].split(":")[1].split("-")
#                     rna_loc = (int(rna_start) + int(rna_end) ) /2
#                     if anchor_chr == rna_chr and abs(int(anchor_loc) - int(rna_loc)) <= 100000:
#                         reported_match.update(["Anchor", "RNASeq"])     
#             # Vs Chm13_Y
#             if not chm13 == "NA":
#                 for c in chm13:
#                     chm13_chr, chm13_pos = c.split(":")
#                     if anchor_chr == chm13_chr and abs(int(anchor_loc) - int(chm13_pos)) <= 1000:
#                         reported_match.update(["Anchor", "Chm13"])
#         #
#         if reported_match:
#             MATCH_dict.append(reported_match)
#     # 10X
#     if (not "10X" in reported_match) and (not x10 == "NA"):
#         reported_match = set()
#         x10_chr, x10_loc = x10.split(":")[:2]
#         # Vs PMC
#         if not pmc == "NA":
#             pmc_fl, pmc_chr, pmc_loc = pmc.split(":")[:3]
#             if x10_chr == pmc_chr and abs(int(x10_loc)*1000 - int(pmc_loc)) <= 100000:
#                 reported_match.update(["10X", "PMC"])
#         # Vs RNASeq
#         if not rna == "NA":
#             for r in rna.split(","):
#                 rna_chr = r.split("|")[0].split(":")[0]
#                 rna_start, rna_end = r.split("|")[0].split(":")[1].split("-")
#                 rna_loc = (int(rna_start) + int(rna_end) ) /2
#                 if x10_chr == rna_chr and abs(int(x10_loc)*1000 - int(rna_loc)) <= 100000:
#                     reported_match.update(["10X", "RNASeq"])
#         # Vs Chm13_Y
#         if not chm13 == "NA":
#             for c in chm13:
#                 chm13_chr, chm13_pos = c.split(":")
#                 if x10_chr == chm13_chr and abs(int(x10_loc)*1000 - int(chm13_pos)) <= 100000:
#                     reported_match.update(["10X", "Chm13"])
#         #
#         if reported_match:
#             MATCH_dict.append(reported_match)
#     # PMC
#     if (not "PMC" in reported_match) and (not pmc == "NA"):
#     #if not pmc == "NA":
#         reported_match = set()
#         pmc_fl, pmc_chr, pmc_loc = pmc.split(":")[:3]
#         # Vs RNASeq
#         if not rna == "NA":
#             for r in rna.split(","):
#                 rna_chr = r.split("|")[0].split(":")[0]
#                 rna_start, rna_end = r.split("|")[0].split(":")[1].split("-")
#                 rna_loc = (int(rna_start) + int(rna_end) ) /2
#                 if rna_chr == pmc_chr and abs(int(rna_loc) - int(pmc_loc)) <= 100000:
#                     reported_match.update(["PMC", "RNASeq"]) 
#         # Vs Chm13_Y
#         if not chm13 == "NA":
#             for c in chm13:
#                 chm13_chr, chm13_pos = c.split(":")
#                 if pmc_chr == chm13_chr and abs(int(pmc_loc) - int(chm13_pos)) <= 1000:
#                     reported_match.update(["PMC", "Chm13"])
#         #
#         if reported_match:
#             MATCH_dict.append(reported_match)
#     # Chm13_Y
#     if (not "Chm13" in reported_match) and (not chm13 == "NA"):
#     #if not chm13 == "NA":
#         reported_match = set()
#         for c in chm13:
#             chm13_chr, chm13_pos = c.split(":")
#             # Vs RNASeq
#             if not rna == "NA":
#                 for r in rna.split(","):
#                     rna_chr = r.split("|")[0].split(":")[0]
#                     rna_start, rna_end = r.split("|")[0].split(":")[1].split("-")
#                     rna_loc = (int(rna_start) + int(rna_end) ) /2
#                     if chm13_chr == rna_chr and abs(int(chm13_pos) - int(rna_loc)) <= 100000:
#                             reported_match.update(["Chm13", "RNASeq"])
#         #
#         if reported_match:
#             MATCH_dict.append(reported_match)

# data = from_memberships(MATCH_dict)

# ## 
# print("Create upsetplot #2")

# fig = plt.figure(figsize = (10, 6)) # width, height
# ax = fig.add_subplot()

# Upset_Intersections = UpSet(data,
#         sort_by = "cardinality", 
#         subset_size = "count", 
#         show_counts = True
#         )

# params = {'font.size': 8}
# with plt.rc_context(params):
#     Upset_Intersections.plot()

# plt.text(-0.8, 0.35, "NRS with matched locations", 
#         fontsize = 10, 
#         transform = ax.transAxes,
#         backgroundcolor = 'white'
#         )


# plt.suptitle("SABE_1172_UNHESMSV")
# plt.savefig("04_SABE_UNHESMSV_Genmap_RNASeq_Chm13Y_Upset_Matches.png")
