import collections
import matplotlib.pyplot as plt
import pandas as pd
import upsetplot
import numpy as np
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import FancyBboxPatch
import re

####
print("Get tissues from samples")
TISSUE_dict = collections.defaultdict(tuple)
DOMAIN_dict = collections.defaultdict(str)

for line in open("07_SABE_1172_UNHESMSV_RNASeq_sample_informations_v2.txt", "r"):
    if line.startswith("#"):
        continue
    Domain, RNASample, Organism, Tissue, CellType, Disease = line.strip("\n").split(";")
    TISSUE_dict[RNASample] = (Tissue, CellType, Disease)
    DOMAIN_dict[RNASample] = Domain

####
print("Total reads per sample | per domain/tissue")
READS_total_dict = collections.defaultdict(int)
READS_DOM_total_dict = collections.defaultdict(int)
all_sampleids, all_domains = set(), set()

for line in open("01_RNASeq_LOG_merged.txt", "r"):
    if line.startswith("#"):
        continue
    sampleid, total_reads = line.split("\t")[:2]
    uniq_mapped = line.split("\t")[4]
    if not uniq_mapped == "0.00%":
        READS_total_dict[sampleid] = int(total_reads)
        all_sampleids.add(sampleid)
        domain = DOMAIN_dict[sampleid]
        all_domains.add(domain)
        READS_DOM_total_dict[domain] += int(total_reads) 

####
print("Reads to NRS per sample")
READS_nrs_dict = collections.defaultdict(int)
READS_DOM_NRS_dict = collections.defaultdict(int)
NRS_dict = collections.defaultdict(set)
NRS_DOM_dict = collections.defaultdict(set)
expressed_nrs = set()

for line in open("04_SABE_1172_UNHESMSV_RNASeq_freq.txt", "r"):
    if line.startswith("#"):
        continue
    nrs, freq, hits = line.strip("\n").split("\t")
    for hit in hits.split(","):
        sampleid, nrs_reads = hit.split(":")
        READS_nrs_dict[sampleid] += int(nrs_reads)
        NRS_dict[sampleid].add(nrs)
        all_sampleids.add(sampleid)
        domain = DOMAIN_dict[sampleid]
        READS_DOM_NRS_dict[domain] += int(nrs_reads) 
        NRS_DOM_dict[domain].add(nrs)
    expressed_nrs.add(nrs)

#### FIGURE FOR PAPER
## Prepare data
# Biotype extraction for expressed NRS
biotype_df = pd.read_csv(
    "06b_Biotypes.txt",
    sep="\t",
    header=None,
    names=["biotype", "count"]
)

# Map the tissues
def map_tissue_to_category(tissue):
    if tissue in {"blood", "bone marrow"}:
        return "Blood / Immune"
    elif tissue == "liver":
        return "Liver"
    elif tissue == "brain":
        return "Brain"
    elif tissue == "skin":
        return "Skin"
    elif tissue == "lung":
        return "Lung"
    else:
        return "Other"
    
NRS_TISSUE_dict = collections.defaultdict(set)

for sampleid, nrs_set in NRS_dict.items():
    tissue_raw = TISSUE_dict[sampleid][0]
    tissue = map_tissue_to_category(tissue_raw)
    for nrs in nrs_set:
        NRS_TISSUE_dict[tissue].add(nrs)

plot_data = []

for tissue, nrs_set in NRS_TISSUE_dict.items():
    plot_data.append({
        "Tissue": tissue,
        "NRS_count": len(nrs_set)
    })

df_plot = pd.DataFrame(plot_data)
df_plot = df_plot.sort_values(by="NRS_count", ascending=False)

##
fig, ax = plt.subplots(figsize=(10, 8))

# Main bar graph
ax.bar(
    df_plot["Tissue"],
    df_plot["NRS_count"],
    color="#6FAFD9",
    edgecolor="black",
    linewidth=0.8
    )

ax.set_xlabel('Tissue', size = 12)
ax.set_ylabel("Number of expressed NRS", size = 12)

ax.tick_params(axis="both", labelsize=10)

ax.set_title(
    "Transcription of NRS across human tissues",
    fontsize=14,
    pad = 10
    )

# Piechart insert
pie_colors = [
    "#FCE5CD",  # Protein-coding (muted salmon)
    "#F4A3A3",  # lncRNA (muted blue)
    "#D9C2E9",  # Pseudogene (muted green)
    "#B6D7A8",  # Small RNA (muted purple)
    "#8FBBD9",  # Immune loci / Other (muted beige)
]
pie_edgecolor = "#626262"

ax_inset = inset_axes(
    ax,
    width="45%",
    height="45%",
    bbox_to_anchor=(0, -0.15, 1, 1),
    bbox_transform=ax.transAxes,
    loc="upper right"
    )

sizes = biotype_df["count"].values
labels = biotype_df["biotype"].values
total = sizes.sum()

wedges, texts, autotexts = ax_inset.pie(
    sizes,
    labels=[
        lbl if (cnt / total) > 0.03 else ""
        for lbl, cnt in zip(labels, sizes)
    ],
    colors=pie_colors[:len(sizes)],
    startangle=90,
    autopct=lambda p: f"{p:.1f}%" if p > 3 else "",
    pctdistance=0.65,
    wedgeprops=dict(edgecolor=pie_edgecolor, linewidth=1),
    textprops=dict(fontsize=10)
)

# Smart label placement for small slices (collision-safe, minimal change)
small_entries = []
for wedge, label, size in zip(wedges, labels, sizes):
    frac = size / total
    if frac < 0.03:
        angle = (wedge.theta2 + wedge.theta1) / 2
        x = np.cos(np.deg2rad(angle))
        y = np.sin(np.deg2rad(angle))

        small_entries.append({
            "label": f"{label} ({frac*100:.1f}%)",
            "x": x,
            "y": y
        })

small_entries.sort(key=lambda d: d["y"])

min_sep = 0.12
for i in range(1, len(small_entries)):
    if small_entries[i]["y"] - small_entries[i-1]["y"] < min_sep:
        small_entries[i]["y"] = small_entries[i-1]["y"] + min_sep

for entry in small_entries:
    x, y = entry["x"], entry["y"]
    ax_inset.annotate(
        entry["label"],
        xy=(x * 0.5, y * 0.6),
        xytext=(x * 1.5, y * 1.35),
        ha="left" if x > 0 else "right",
        va="center",
        fontsize=10,
        arrowprops=dict(
            arrowstyle="-",
            color="#727272",
            lw=0.8
        )
    )

ax_inset.set_title(
    "Genomic context of expressed NRS",
    fontsize=12,
    pad = 50
    )


bbox = FancyBboxPatch(
    (-0.1, 0.05), 1.2, 1.25,
    boxstyle="round,pad=0.03",
    linewidth=1.0,
    edgecolor="#D7D7D7F8",
    facecolor="white",
    transform=ax_inset.transAxes,
    zorder=0,
    clip_on = False
)
ax_inset.add_patch(bbox)

# Hide default spines
for spine in ax_inset.spines.values():
    spine.set_visible(False)


# Plot parameters
#fig.tight_layout()
fig.savefig("Fig4B_NRS_expression_across_tissues.png", dpi=300)
plt.close(fig)

####
# print("Print to file - per Domain")

# outfile = open("08_Reads_per_Domain.txt", "w+")
# outfile.write("#Domain\tTotal_Reads\tNRS_Reads\tNRS_perc\tpercentage\tNRS_list\n")

# x_domain_list, y_percentage_list = list(), list()
# plot_data = []

# for domain in all_domains:
#     total_reads = READS_DOM_total_dict[domain]
#     NRS_reads = READS_DOM_NRS_dict[domain]
#     percentage = round((NRS_reads * 100)/ total_reads, 2)
#     nrs_list = NRS_DOM_dict[domain]
#     #
#     plot_data.append({"Domain": domain,
#                     "percentage": percentage})
#     # x_domain_list.append(domain)
#     # y_percentage_list.append(percentage)
#     #
#     outfile.write("\t".join([domain, str(total_reads), str(NRS_reads), str(percentage), str(len(nrs_list)), ",".join(nrs_list)]) + "\n")

# outfile.close()

# df = pd.DataFrame(plot_data)
# df = df.sort_values(by = "percentage", axis = 0, ascending = True)

# # Merge the tissue types into tissue categories
# outfile = open("08_Reads_per_Sample.txt", "w+")
# outfile.write("#NRS\tTotal_Reads\tNRS_Reads\tTotal_NRS\tNRS_list\n")

# list_of_selected_tissues = ["blood", "heart", "liver", "brain",  "skin", "intestines"]

# data = []
# data2 = collections.defaultdict(set)

# for sampleid in all_sampleids:
#     total_reads = READS_total_dict[sampleid]
#     nrs_reads = READS_nrs_dict[sampleid]
#     nrs_list = NRS_dict[sampleid]
#     tissue = TISSUE_dict[sampleid][0]
#     disease = TISSUE_dict[sampleid][2]
#     cell_type = TISSUE_dict[sampleid][1]
#     ##
#     outfile.write("\t".join([sampleid, str(total_reads), str(nrs_reads), str(len(nrs_list)), ",".join(nrs_list)]) + "\n")
#     ##
#     data.append({"SampleID": sampleid, 
#                 "Total_Reads_%": (total_reads-nrs_reads)*100/total_reads, 
#                 "NRS_reads_%": nrs_reads*100/total_reads, 
#                 "NRS_count": len(nrs_list),
#                 "Tissue": tissue,
#                 "Disease": disease,
#                 "Cell_Type": cell_type
#                 })
#     ##
#     if disease == "normal" and tissue in list_of_selected_tissues:
#         for nrs in nrs_list:
#             data2[tissue].add(nrs)

# outfile.close()


#### UPSETPLOT disabled May-2026 due to packages version incompatibility
#####

## Prepare for graphs
# df = pd.DataFrame(data)

# df = df.sort_values(by = "NRS_reads_%", axis = 0, ascending = True) 
# print("Top 10 tissues sorted by NRS reads %:")
# print(df.head(10))


# df = df.sort_values(by = "NRS_count", axis = 0, ascending = False) 
# print("Top 10 tissues sorted by NRS count:")

### Upset plot for NRS in tissues
#df2 = pd.DataFrame(data2)
# d2 = upsetplot.from_contents(data2)

# fig, ax = plt.subplots()

# data2 = upsetplot.plot(d2,
#         sort_by = "cardinality",
#         #subset_size = "sum",
#         #sum_over = "length",
#         show_counts = '%.1f'
#         )

# #plt.suptitle("SABE_1172_UNHESMSV")
# plt.savefig("08_NRS_across_tissue_upsetplot.png")

#### Potential supplementary graphs
# ## Reads aligned to NRS histogram

# fig = plt.figure() #figsize = (10, 8))
# ax = fig.add_subplot()

# ax.hist(df["NRS_reads_%"], bins = 30, log = True)
# ax.set_xlabel("Reads aligned to NRS (%) per sample")

# fig.tight_layout()
# plt.savefig("08_Reads_per_sampleid_hist.png")

# ##
# fig = plt.figure() #figsize = (10, 8))
# ax = fig.add_subplot()


# ax.hist(df["NRS_count"], bins = 40, log = True)
# ax.set_xlabel("NRS count per sample")

# fig.tight_layout()
# plt.savefig("08_NRS_per_sampleid_hist.png")

## Reads aligned to NRS (%) histogram

# fig = plt.figure() #figsize = (10, 8))
# ax = fig.add_subplot()

# ax.hist(df["percentage"], 
#         bins = 30,
#         weights = np.ones(len(df["percentage"]))/len(df["percentage"]))

# ax.set_xlabel("NRS per domain (%)")
# ax.set_ylabel("NRS fraction")

# fig.tight_layout()
# plt.savefig("08_Reads_per_domain_hist.png")
