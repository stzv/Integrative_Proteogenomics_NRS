
import pandas as pd
import matplotlib.pyplot as plt
from math import log
import numpy as np
import collections
import seaborn as sns

import sys

####


####
NRS_dict = collections.defaultdict(dict)
Tissues_dict = collections.defaultdict(set)
Frequency_dict = collections.defaultdict(dict)

for line in open("07_SABE_1172_UNHESMSV_NRS_RNA_tissue_table.txt", "r"):
    if line.startswith("#"):
        continue
    ##
    NRS, Popul_Freq, RNA_Freq, Tissues, Exons = line.strip("\n").split("\t")
    nrs_split = NRS.split("_")
    NRS_len = int(nrs_split[-1]) - int(nrs_split[-2]) + 1
    ##
    if RNA_Freq == "0":
        continue
    ##
    for Tiss in Tissues.split(","):
        if not Tiss.startswith(("ERR", "SRR")):
            continue
        rna_sample, cpm, tissue = Tiss.split(":")
        Tissues_dict[tissue].add((NRS, NRS_len, Popul_Freq, rna_sample, cpm))

## Perform manual magic of selecting which tissue belongs to which category
categories = dict()
for line in open("07a_categories_tissue_manual.txt", "r"):
    category = line.split(";")[0]
    tissue_names = line.strip("\n").split(";")[1:]
    for tissue in tissue_names:
        categories[tissue] = category

##
Tissue_AvgCPM_dict = collections.defaultdict(float)
categories_samples = collections.defaultdict(set)
categories_nrs_count = collections.defaultdict(set)
categories_nrs_freq = collections.defaultdict(list)

for Tissue, entries in Tissues_dict.items():
    cpms = list()
    for entry in entries:
        nrs, bp, Popul_Freq, RNA, CPM = entry
        # Only take samples where the expression is above 0.2
        if float(CPM) < 0.2:
            continue
        cpms.append(float(CPM))
        # Get the numbers for tissues that are part of selected categories
        if Tissue in categories.keys():
            categories_samples[categories[Tissue]].add(RNA)
            categories_nrs_count[categories[Tissue]].add(nrs)
            categories_nrs_freq[categories[Tissue]].append(float(Popul_Freq))
    ##
    average_cpm = sum(cpms) / len(cpms)
    Tissue_AvgCPM_dict[Tissue] = round(average_cpm, 2)
    ##
    
### Print out the average CPM into file
categories_cpm = collections.defaultdict(list)

with open("07a_all_tissues_average_cpm.txt", "w") as outfile:
    outfile.write("#Tissue\tMeanCPM\tCategory")
    for Tissue, cpm in Tissue_AvgCPM_dict.items():
        if Tissue in categories.keys():
            category = categories[Tissue]
            #if float(cpm) < 200:
            categories_cpm[category].append(log(float(cpm)))
        else:
            category = ""
        outfile.write(f"{Tissue}\t{cpm}\t{category}\n")

# CPM distribution
data = pd.concat([pd.DataFrame({k: v}) for k, v in categories_cpm.items() if len(v) >= 5], axis = 1).melt().dropna()
# RNASeq samples count
data2 = pd.DataFrame({k: len(v)} for k, v in categories_samples.items() if k in list(data["variable"])).melt().dropna()
# Expressed NRS count
data3 = pd.DataFrame({k: (len(v))} for k, v in categories_nrs_count.items() if k in list(data["variable"])).melt().dropna()
# Population frequency of expressed NRS
data4 = pd.concat([pd.DataFrame({k: v}) for k, v in categories_nrs_freq.items() if k in list(data["variable"])], axis = 1).melt().dropna()

####
fig = plt.figure(figsize = (10, 8))

#Plot RNASeq samples count per category
ax1 = fig.add_subplot(111)
sns.swarmplot(data = data, x = 'variable', y = 'value', size = 4)#, color = "black")


# Plot distribution of CPM per category
ax2 = ax1.twinx()
ax2 = sns.barplot(data = data3, x = "variable", y = "value",
            linewidth = 1,
            #color = "none"
            alpha = 0.3
            )


# Add the labels to the barplot
rects = ax2.patches
labels = list(data3["value"])

for rect, label in zip(rects, labels):
    height = rect.get_height()
    ax2.text(
        rect.get_x() + rect.get_width() / 2, 
        height + 5, 
        f"{label:,.10g}", 
        ha = "center", va = "bottom",
        size = "small"
    )


# Push the barplot to the background
ax1.set_zorder(ax2.get_zorder()+1)
ax1.patch.set_visible(False) # type: ignore

# Adjust the visual of the axis
x_ticks =  [item.get_text() for item in ax1.get_xticklabels()]
x_ticks_new = list()

for x in x_ticks:
    rna = len(categories_samples[x])
    x_new = f"{x}\n({rna:,.5g})"
    x_ticks_new.append(x_new)


ax1.set_xticklabels(x_ticks_new, rotation = 45, ha = "center")

##

ax2.ticklabel_format(axis = "y", style = "scientific", scilimits = (4, 4))

##

ax1.set_ylabel("Average log CPM per tissue")
ax1.set(xlabel = "Tissue categories (Samples count)")

ax2.set_ylabel("Expressed NRS count")

#sns.lineplot(data = data2, x = "variable", y = "value", label = "RNASeq sample count", alpha = 0.2)
# ax2.set_ylabel("RNASeq sample count per category")


# # Plot NRS count per category
# ax3 = fig.add_subplot(212)
# ax3 = sns.barplot(data = data3, x = "variable", y = "value", label = "NRS count")

# # for i in ax3.containers:
# #     ax3.bar_label(i,)

# ax3.xaxis.tick_top()
# y_ticks =  ax3.get_yticks()
# ax3.set_yticklabels([int(abs(tick)) for tick in y_ticks])

# ax3.set_ylabel("Expressed NRS")
# ax3.set(xlabel = None)
# ax3.set(xticklabels = [])


fig.tight_layout()
plt.savefig("07a_SABE_1172_UNHESMSV_Expression_Picture.png")


