import pandas as pd
import sys
import collections
from matplotlib import pyplot as plt
from upsetplot import from_memberships
from upsetplot import UpSet
import upsetplot
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import matplotlib.colors as mcolors


#####
print("Compare locations #1")
RNASeq_dict = collections.defaultdict(dict)

for line in open("06_SABE_1172_UNHESMSV_NRS_best_CJ.txt", "r"):
    nrs, locs = line.strip("\n").split("\t")
    nrs_short = "_".join(nrs.split("_")[:2])
    if locs.startswith("HLA"):
        continue
    RNASeq_dict[nrs_short] = locs

###### Upsetplot - only whether they have location or not

NRS_dict = collections.defaultdict(int)
LOC_dict = collections.defaultdict(set)

outfile = open("02_SABE_UNHESMSV_Merged_Locations.txt", "w+")
outfile.write("#NRS\tLength_bp\tPMC\tX10\tAnchors\tRNASeq\n")


for line in open("01b_SABE_1172_UNHESMSV_locations2.txt", "r"):
    nrs, bp, pmc, x10, anchors = line.strip("\n").replace("chr", "").split("\t")
    rna = "NA"
    if nrs in RNASeq_dict.keys():
        rna = RNASeq_dict[nrs]
    NRS_dict[nrs] = int(bp)/1000000
    ## Print to file
    if not x10 == "NA":
        x10_chr, x10_loc, x10_reads = x10.split(":")
        X10 = ":".join([x10_chr, str(int(x10_loc)*1000), x10_reads])
        outline = "\t".join([nrs, bp, pmc, X10, anchors, rna])
    else:
        outline = "\t".join([nrs, bp, pmc, x10, anchors, rna])
    outfile.write(outline + "\n")
    ##
    if not pmc == "NA":
        LOC_dict["PMC"].add(nrs)
    if not anchors == "NA":
        LOC_dict["Anchors"].add(nrs)
    if not rna == "NA":
        LOC_dict["RNASeq"].add(nrs)
    if not x10 == "NA":
        LOC_dict["10X"].add(nrs)
    ## No location at all
    if all(x  == "NA" for x in (pmc, x10, anchors, rna)):
        LOC_dict["Unloc"].add(nrs)

### preparation for upsetplot

nrs_df = pd.DataFrame(list(NRS_dict.items()))
nrs_df.columns = ["nrs", "length"]

anch_df = pd.DataFrame(list(LOC_dict["Anchors"]))
anch_df.columns = ["nrs"]
anch_df["length"] = anch_df['nrs'].map(NRS_dict)

pmc_df = pd.DataFrame(list(LOC_dict["PMC"]))
pmc_df.columns = ["nrs"]
pmc_df["length"] = pmc_df['nrs'].map(NRS_dict)

x10_df = pd.DataFrame(list(LOC_dict["10X"]))
x10_df.columns = ["nrs"]
x10_df["length"] = x10_df['nrs'].map(NRS_dict)

rna_df = pd.DataFrame(list(LOC_dict["RNASeq"]))
rna_df.columns = ["nrs"]
rna_df["length"] = rna_df['nrs'].map(NRS_dict)

noloc_df = pd.DataFrame(list(LOC_dict["Unloc"]))
noloc_df.columns = ["nrs"]
noloc_df["length"] = rna_df['nrs'].map(NRS_dict)

### !!! Generating upsetplot disabled on 07-May-2026 due to incompatibility between python 3.13, pandas > 3.0 and upsetplot.
# # Draw the uspetplot
# print("Create upsetplot #1")

# content2 = {"Total NRS": [nrs_df["nrs"], nrs_df["length"]],
#         "Anchors": [anch_df["nrs"], anch_df["length"]],
#         "PMC": [pmc_df["nrs"], pmc_df["length"]],
#         "10X": [x10_df["nrs"], x10_df["length"]],
#         "RNASeq": [rna_df["nrs"], rna_df["length"]]     
#         }

# data2 = upsetplot.from_contents(content2)


# #
# plt.rcParams["scatter.edgecolors"] = "black"
# fig = plt.figure(figsize = (7, 6)) # width, height
# ax = fig.add_subplot()

# Upset_Intersections = UpSet(data2,
#         sort_by = "cardinality", subset_size = "count", 
#         show_counts = True,
#         )

# params = {'font.size': 8}
# with plt.rc_context(params):
#     Upset_Intersections.plot()

# plt.text(-0.1, 0.35, "Mapped NRS count", 
#         fontsize = 10, 
#         transform = ax.transAxes,
#         backgroundcolor = 'white'
#         )

# plt.suptitle("SABE_1172_UNHESMSV")
# plt.savefig("02_SABE_UNHESMSV_Genmap_RNASeq_Upset_LocationsCount.png")


# -----------------------------
# Contribution of mapping methods
# (Figure 3B)
# -----------------------------

# Total number of NRS
total_nrs = len(NRS_dict)

# Counts per method
pmc_set     = (LOC_dict["PMC"])
rnaseq_set  = (LOC_dict["RNASeq"])
anchor_set  = (LOC_dict["Anchors"])
x10_set     = (LOC_dict["10X"])

# Cumulative unions 
anchor_only = anchor_set
rnaseq_new = rnaseq_set - anchor_set
x10_new = x10_set - (anchor_set | rnaseq_set)
pmc_new = pmc_set - (anchor_set | rnaseq_set | x10_set)
placed_all = anchor_set | rnaseq_set | x10_set | pmc_set
unmapped = set(NRS_dict.keys()) - placed_all

# Labels for legend
n_anchor = len(anchor_only)
n_rnaseq = len(rnaseq_new)
n_x10    = len(x10_new)
n_pmc    = len(pmc_new)
n_unmap  = len(unmapped)

# Fractions
frac_anchor = len(anchor_only) / total_nrs
frac_rnaseq = len(rnaseq_new) / total_nrs
frac_x10    = len(x10_new) / total_nrs
frac_pmc    = len(pmc_new) / total_nrs
frac_unmap  = len(unmapped) / total_nrs

labels = [
    "Anchoring",
    "+ RNA‑Seq",
    "+ 10x linked‑reads",
    "+ PMC",
    "All methods"
]
x = np.arange(len(labels))

# Colors
col_anchor = "#fd8d3c"
col_rnaseq = "#74c476"
col_x10    = "#6baed6"
col_pmc    = "#ac5eb4"
col_unmap  = "#a6a6a6"

# Plot

fig, ax = plt.subplots(figsize=(7.5, 4.8))

# Color for placed fraction
placed_color = "#6baed6"   # consistent blue
unmapped_color = "#bdbdbd" # neutral gray

# Plot cumulative placed fraction
text_offset = 0.015
# Bar 1: Anchoring
ax.bar(x[0], frac_anchor, color=col_anchor, edgecolor="black")
ax.text(x[0], frac_anchor + text_offset, f"{frac_anchor:.1}\n({n_anchor:,})",
    ha="center", va="bottom", fontsize=9)

# Bar 2: Anchoring + RNA-Seq
ax.bar(x[1], frac_anchor, color=col_anchor, edgecolor="black")
ax.bar(x[1], frac_rnaseq, bottom=frac_anchor,
        color=col_rnaseq, edgecolor="black")
ax.text(x[1], frac_anchor + frac_rnaseq + text_offset, f"+{frac_rnaseq:.1}\n({n_rnaseq:,})",
    ha="center", va="bottom", fontsize=9)

# Bar 3: + 10x
ax.bar(x[2], frac_anchor, color=col_anchor, edgecolor="black")
ax.bar(x[2], frac_rnaseq, bottom=frac_anchor,
        color=col_rnaseq, edgecolor="black")
ax.bar(x[2], frac_x10, bottom=frac_anchor + frac_rnaseq,
        color=col_x10, edgecolor="black")
ax.text(x[2], frac_anchor + frac_rnaseq + frac_x10 + text_offset, f"+{frac_x10:.1}\n({n_x10:,})",
    ha="center", va="bottom", fontsize=9)

# Bar 4: + PMC
ax.bar(x[3], frac_anchor, color=col_anchor, edgecolor="black")
ax.bar(x[3], frac_rnaseq, bottom=frac_anchor,
        color=col_rnaseq, edgecolor="black")
ax.bar(x[3], frac_x10, bottom=frac_anchor + frac_rnaseq,
        color=col_x10, edgecolor="black")
ax.bar(x[3], frac_pmc,
        bottom=frac_anchor + frac_rnaseq + frac_x10,
        color=col_pmc, edgecolor="black")
ax.text(x[3], frac_anchor + frac_rnaseq + frac_x10 + frac_pmc + text_offset, f"+{frac_pmc:.1}\n({n_pmc:,})",
    ha="center", va="bottom", fontsize=9)

# Bar 5: All + unmapped
ax.bar(x[4], frac_anchor, color=col_anchor, edgecolor="black", label=f"Anchoring")
ax.bar(x[4], frac_rnaseq, bottom=frac_anchor,
        color=col_rnaseq, edgecolor="black", label=f"RNA‑Seq chimeric")
ax.bar(x[4], frac_x10, bottom=frac_anchor + frac_rnaseq,
        color=col_x10, edgecolor="black", label=f"10x linked‑reads")
ax.bar(x[4], frac_pmc,
        bottom=frac_anchor + frac_rnaseq + frac_x10,
        color=col_pmc, edgecolor="black", label=f"Partial mapping (PMC)")
ax.bar(x[4], frac_unmap,
        bottom=frac_anchor + frac_rnaseq + frac_x10 + frac_pmc,
        color=col_unmap, edgecolor="black", label=f"Unmapped")
ax.text(x[4], frac_anchor + frac_rnaseq + frac_x10 + frac_pmc + (frac_unmap / 2), f"{frac_unmap:.1}\n({n_unmap:,})",
    ha="center", va="bottom", fontsize=9)

# Plot settings
ax.set_ylabel("Fraction of NRS", fontsize=11)
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=10)
ax.set_ylim(0, 1.0)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.legend(frameon=False, loc="upper left")

plt.tight_layout()
plt.savefig("Fig3_NRS_Placement_Method_Contribution.png", dpi=300)
plt.close()

########################
## In case the rest doesn't need to run:
import sys
sys.exit()
########################


###### Upsetplot - match the locations
print("Compare locations #2")

NRS_dict = collections.defaultdict(int)
LOC_dict = collections.defaultdict(set)

MATCH_dict = list()
NRS_matches_dict = collections.defaultdict(list)

output = open("02_NRS_GENMAP.tsv", "w")
output.write("#NRS\tLength_bp\tGENMAP\tConsens_Loc\tGENMAP_Method\tMatches\tPMC\t10X\tAnchors\tRNASeq\n")

counter = 0


for line in open("01b_SABE_1172_UNHESMSV_locations2.txt", "r"):
    ## Prepare variables
    nrs, bp, pmc, x10, anchors = line.strip("\n").replace("chr", "").split("\t")
    ##
    rna = "NA"
    if nrs in RNASeq_dict.keys():
        rna = RNASeq_dict[nrs]
    NRS_dict[nrs] = int(bp)/1000000
    ## Match the locations
    matched_anchor = ""
    matched_rna = ""
    # 10X
    if not x10 == "NA": 
        reported_match = set()
        x10_chr, x10_loc = x10.split(":")[:2]
        # Vs Anchor
        if not anchors == "NA":
            for anchor in anchors.split(","):
                anchor_chr, anchor_loc = anchor.split(":")[:2]  
                if anchor_chr == x10_chr and abs(int(anchor_loc) - (int(x10_loc)*1000)) <= 50000:
                    reported_match.update(["10X", "Anchor"])
                    matched_anchor = anchor
        # Vs PMC
        if not pmc == "NA":
            pmc_fl, pmc_chr, pmc_loc = pmc.split(":")[:3]
            if x10_chr == pmc_chr and abs(int(x10_loc)*1000 - int(pmc_loc)) <= 50000:
                reported_match.update(["10X", "PMC"])
        # Vs RNASeq
        if not rna == "NA":
            for r in rna.split(","):
                rna_chr = r.split("|")[0].split(":")[0]
                rna_start, rna_end = r.split("|")[0].split(":")[1].split("-")
                rna_loc = (int(rna_start) + int(rna_end) ) /2
                if x10_chr == rna_chr and abs(int(x10_loc)*1000 - int(rna_loc)) <= 50000:
                    reported_match.update(["10X", "RNASeq"])
                    matched_rna = r
        #
        if reported_match:
            MATCH_dict.append(reported_match)
            NRS_matches_dict[nrs].append(reported_match)
    # Anchor
    if (not "Anchor" in [i for item in NRS_matches_dict[nrs] for i in item]) and (not "NA" in anchors):
        reported_match = set()
        for anchor in anchors.split(","):
            anchor_chr, anchor_loc = anchor.split(":")[:2]
            # Vs PMC
            if not pmc == "NA":
                pmc_fl, pmc_chr, pmc_loc = pmc.split(":")[:3]
                if anchor_chr == pmc_chr and abs(int(anchor_loc) - int(pmc_loc)) <= 1000:
                    reported_match.update(["Anchor", "PMC"])
                    matched_anchor = anchor
            # Vs RNASeq
            if not rna == "NA":
                for r in rna.split(","):
                    rna_chr = r.split("|")[0].split(":")[0]
                    rna_start, rna_end = r.split("|")[0].split(":")[1].split("-")
                    rna_loc = (int(rna_start) + int(rna_end) ) / 2
                    if anchor_chr == rna_chr and abs(int(anchor_loc) - int(rna_loc)) <= 50000:
                        reported_match.update(["Anchor", "RNASeq"])
                        matched_rna = r
                        matched_anchor = anchor
        #
        if reported_match:
            MATCH_dict.append(reported_match)
            NRS_matches_dict[nrs].append(reported_match)
    # PMC
    if (not "PMC" in [i for item in NRS_matches_dict[nrs] for i in item]) and (not pmc == "NA"):
        reported_match = set()
        pmc_fl, pmc_chr, pmc_loc = pmc.split(":")[:3]
        # Vs RNASeq
        if not rna == "NA":
            for r in rna.split(","):
                rna_chr = r.split("|")[0].split(":")[0]
                rna_start, rna_end = r.split("|")[0].split(":")[1].split("-")
                rna_loc = (int(rna_start) + int(rna_end) ) /2
                if rna_chr == pmc_chr and abs(int(rna_loc) - int(pmc_loc)) <= 50000:
                    reported_match.update(["PMC", "RNASeq"])
                    matched_rna = r
        #
        if reported_match:
            MATCH_dict.append(reported_match)
            NRS_matches_dict[nrs].append(reported_match)
    ## Print the consensus sequence
    consensus = "NA"
    method = "NA"
    matches =  [":".join(match) for match in NRS_matches_dict[nrs]]
    if len(matches) == 1: # If only 1 matching set
        if "10X" in str(matches): #10X most reliable, but imprecise
            consensus = x10
            method = "10X Linked Reads"
            if "Anchor" in str(matches):
                consensus = matched_anchor
                method = "Anchoring"
            elif "RNASeq" in str(matches):
                consensus = matched_rna
                method = "Chimeric Reads"
        elif "Anchor" in str(matches): # No 10X but anchor in match
            consensus = matched_anchor
            method = "Anchoring"
        elif str(matches) == "PMC:RNASeq": # No 10 or Anchor, but PMC match RNASeq
            consensus = matched_rna
            method = "Chimeric Reads"
    else: # If no match or two different pairs of match
        if not x10 == "NA": # take 10X if any
            consensus = x10
            method = "10X Linked Reads"
        elif not anchors == "NA": #Anchor if only 1 anchor
            if not len(anchors.split(",")) > 1:
                consensus = anchors
                method = "Anchoring"
    # 
    matches = ",".join(matches)
    if not matches:
        matches = "NA"
    # Rewrite consensus format 
    if not (consensus.startswith("M") or consensus.startswith("hs38d1") or "v" in consensus): #if not MT or decoy
        if method == "Anchoring":
            chrom, pos, reads, direction, strand = consensus.split(":")
            start = int(pos) - 500 + 1 # 1-based
            end = int(pos) + 500 
            GENMAP = f"chr{chrom}:{start}:{end}:{strand}"
        elif method == "Chimeric Reads":
            chrom = consensus.split(":")[0]
            start, end = consensus.split(":")[1].split("|")[0].split("-")
            strand = "*"
            GENMAP = f"chr{chrom}:{start}:{end}:{strand}"
        elif method == "10X Linked Reads": # Precision <=> 100kb
            chrom, pos, read_info = consensus.split(":")
            start = (int(pos)*1000) - 50000
            end = (int(pos)*1000) + 50000
            strand = "*"
            GENMAP = f"chr{chrom}:{start}:{end}:{strand}"
        elif method == "NA":
            GENMAP = "NA"
        else:
            print("ERROR")
            print(nrs, consensus)
            sys.exit()
    else:
        GENMAP = consensus
    #
    if consensus.startswith("chrEBV"):
        consensus = consensus.replace("chrEBV", "EBV")
    outline = "\t".join([nrs, bp, GENMAP, consensus, method, matches, pmc, x10, anchors, rna])
    output.write(outline + "\n")
    if not consensus == "NA":
        counter += 1

print("Total mapped", counter)

### !!! Generating upsetplot disabled on 07-May-2026 due to incompatibility between python 3.13, pandas > 3.0 and upsetplot.
# ## 
# print("Create upsetplot #2")
# data = from_memberships(MATCH_dict)

# fig = plt.figure(figsize = (10, 6)) # width, height
# ax = fig.add_subplot()

# Upset_Intersections = UpSet(data,
#                             sort_by = "cardinality", 
#                             subset_size = "count", 
#                             show_counts = True
#                             )

# params = {'font.size': 8}
# with plt.rc_context(params):
#     Upset_Intersections.plot()

# plt.text(-0.12, 0.25, "NRS with matched locations", 
#         fontsize = 10, 
#         transform = ax.transAxes,
#         backgroundcolor = 'white'
#         )


# plt.suptitle("SABE_1172_UNHESMSV")
# plt.savefig("02_SABE_UNHESMSV_Genmap_RNASeq_Upset_Matches.png")

#### Contribution plot for publication
