import sys

##
NRS_dict, X10_dict, ANCHORS_dict, PMC_dict = dict(), dict(), dict(), dict()

##
print("Load in NRS")
for line in open("SABE_1172_UNHESMSV_genotyping/SABE1172_UNHESMSV_NRS_dark_freeze_final.fa", "r"):
    if line.startswith(">"):
        ctg = line.split("_")[0].replace(">", "") + "_" + line.split("_")[1]
        length = int(line.split("_")[3].split(" ")[0]) - int(line.split("_")[2]) + 1
        NRS_dict[ctg] = length

#
print("Load in 10X")
for line in open("02_SABE_UNHESMSV_1172_10Xlocs_readsbc_sum.txt", "r"):
    ctg, top10x, hits = line.split("\t")
    ctg = ctg.split("_")[0] + "_" + ctg.split("_")[1]
    samples = top10x.split("|")[3]
    barcodes = int(top10x.split("|")[2]) # Extract position's barcode sum
    xreads = int(top10x.split("|")[1])
    top10x = top10x.split("|")[0].replace("chr", "").replace("Un_", "")
    X10_dict[ctg] = {"top10x": top10x, "bc": int(barcodes), "reads": int(xreads), "samples": samples}

#
print("Load in anchors")
for line in open("03_SABE_UNHESMSV_1172_anchors_reads_count.txt", "r"):
    ctg, anchors = line.strip("\n").split("\t")
    anchors = anchors.strip("[]").replace("'", "").replace(" ", "").split(",")
    # Extract top anchores
    anchor_scores = list()
    for ansc in anchors:
        sc = ansc.split("|")[0]
        anchor_scores.append(sc)
    max_anchor_score = max(anchor_scores)
    anchs = list()
    for ansc in anchors:
        if ansc.startswith(f"{max_anchor_score}|"):
            anchs.append(ansc)
    ANCHORS_dict[ctg] = anchs

#
print("Load in PMC")
for line in open("SABE_1172_UNHESMSV_genotyping/SABE_1172_UNHESMSV_PMC_locs.txt", "r"):
    ctg, position = line.strip("\n").split("\t")
    ctg = ctg.split("_")[0] + "_" + ctg.split("_")[1]
    chrom, pos = position.split(":")[1:3]
    PMC_dict[ctg] = {"chrom": chrom, "pos": pos}
#
#
print("\nCombine locations")
counter_both_locations, counter_match, counter_notmatch_singletons = 0, 0, 0
supported_pmc_10x, supported_pmc_anchor = 0, 0
not_matched_bc, not_matched_reads, matched_bc, matched_reads, not_matched_10xreads,not_matched_NRS = list(), list(), list(), list(), list(), list()
supported_pmc = list()
not_matched_compare = dict()

match_out = open("05_SABE_1172_UNHESMSV_10x_anchor_match.txt", "w+")
match_out.write("nrs\tPMC\t10x_loc\t10xbc\treads\tanchors_matched\tall_anchors\n")

not_matched_out = open("05_SABE_1172_UNHESMSV_10x_anchor_notmatch.txt", "w+")
not_matched_out.write("NRS\tPMC\tSupports\t10XLoc\tsamples\t10Xbarcode\t10Xreads\tmaxanchorreads\tanchors\n")
    
#only_anchor = open("SABE_1172_UNHESMSV_onlyanchor.txt", "w+")
#only_10x = open("SABE_1172_UNHESMSV_only10X.txt", "w+")

for ctg in NRS_dict:
    pmc_chrom, pmc_pos, pmc_loc, supports = "NA", "NA", "NA", "NA"
    # If ctg doesn't have both anchor and 10X, skip
    if (not ctg in X10_dict.keys()) or (not ctg in ANCHORS_dict.keys()):
        #if ctg in X10_dict.keys():
        #    only_10x.
        continue
    # If have both locations
    counter_both_locations += 1
    x10 = X10_dict[ctg]["top10x"]
    x10_chr, x10_pos = x10.split(":")[:2]
    x10_pos = int(x10_pos) * 1000
    anchor = ANCHORS_dict[ctg]
    # Compare locations
    MATCH = "NO"
    matched_anchors, reads_list, reads_list_matched, anchor_chroms = list(), list(), list(), list()
    for a in anchor: # for every anchor hit for the NRS
        reads = a.split("|")[0]
        hits = a.split("|")[1:]
        for h in hits: # for every anchor within the cluster
            chrom = h.split(":")[0]
            #for PMC check
            if chrom not in anchor_chroms:
                anchor_chroms.append(chrom)
            #
            pos = h.split(":")[1]
            direction = h.split(":")[2]
            strand = h.split(":")[3]
            # For each anchor, compare with 10X
            if (x10_chr == chrom) and (abs(x10_pos - int(pos)) < 100000): #if chrom match and anchor <> 10x pos is less than 200kb (100 each side)
                MATCH = "YES"
                matched_anchors.append(f"chr{chrom}:{pos}:{direction}:{strand}")
                reads_list_matched.append(int(reads))
            else: # if there is no match, save the reads count of the anchors
                reads_list.append(int(reads))
    # Check if locations have a match
    if MATCH == "YES":
        counter_match += 1
        # check with PMC location
        if ctg in PMC_dict.keys():
            pmc_chrom = PMC_dict[ctg]["chrom"]
            pmc_pos = int(PMC_dict[ctg]["pos"])
            pmc_loc = pmc_chrom + ":" + str(pmc_pos)
        # save the match
        match_out.write(f"{ctg}\t{pmc_loc}\t{x10}\t{X10_dict[ctg]['bc']}\t{max(reads_list_matched)}\t{','.join(matched_anchors)}\t{anchor}\n")
        # for graph
        matched_bc.append(X10_dict[ctg]["bc"])
        matched_reads.append(max(reads_list_matched))
        continue # Continue with next NRS
    elif MATCH == "NO":
        # If 10X and anchor did not match for NRS, compare the read count of anchors with barcodes of 10X
        # only singleton anchors count
        if int(max(reads_list)) < 2:
            counter_notmatch_singletons += 1
        # for graph
        not_matched_reads.append(max(reads_list))
        not_matched_bc.append(int(X10_dict[ctg]["bc"]))
        not_matched_10xreads.append(int(X10_dict[ctg]["reads"]))
        not_matched_NRS.append(ctg)
        # Check if 10X or anchor supported by PMC if anchor does not match 10X
        if ctg in PMC_dict.keys():
            supports = "N/A"
            pmc_chrom = PMC_dict[ctg]["chrom"]
            pmc_pos = int(PMC_dict[ctg]["pos"])
            pmc_loc = pmc_chrom + ":" + str(pmc_pos)
            if (pmc_chrom == x10_chr) and (abs(pmc_pos - x10_pos) < 100000):
                supported_pmc_10x += 1
                supported_pmc.append(ctg)
                supports = "10X"
            elif pmc_chrom in anchor_chroms:
                for an in anchor:
                    hit = an.split("|")[1:]
                    matched2 = "NO"
                    for hi in hit:
                        hi_chrom = hi.split(":")[0]
                        hi_pos = hi.split(":")[1]
                        if (pmc_chrom == hi_chrom) and (abs(pmc_pos - int(hi_pos)) < 1000):
                            matched2 = "YES"
            if matched2 == "YES":
                supported_pmc_anchor+= 1
                supported_pmc.append(ctg)
                supports = "ANCH"
                #
                # print out only if anchor reads > 5
                #if max(reads_list) > 4:
            not_matched_out.write(f"{ctg}\t{pmc_loc}\t{supports}\t{x10}\t{X10_dict[ctg]['samples']}\t{X10_dict[ctg]['bc']}\t{X10_dict[ctg]['reads']}\t{max(reads_list)}\t{anchor}\n")
    #

match_out.close()
not_matched_out.close()

#print("\tNot matched, only singleton anchors", counter_notmatch_singletons)
#print("\tMin & max unmatched anchor reads", min(not_matched_reads), max(not_matched_reads))
#print("\tMin & max unmatched barcodes", min(not_matched_bc), max(not_matched_bc))
#print("\tMin & max unmatched 10X reads", min(not_matched_10xreads), max(not_matched_10xreads))
print("\n\tNRS mapped", counter_match)
print("\tNRS not mapped", len(not_matched_NRS))
print("\t\tSupported by PMC:")
print(" \t\t10X", supported_pmc_10x)
print(" \t\tAnchors", supported_pmc_anchor)

sys.exit()

####
print("Print graphs")
import matplotlib.pyplot as plt
import pandas as pd

df_nomatch = pd.DataFrame(data = {'nrs': not_matched_NRS, 'barcode': not_matched_bc, 'x10reads': not_matched_10xreads, 'reads': not_matched_reads}) 

print(df_nomatch.head(5))

## Barcodes
plt.scatter(df_nomatch["barcode"], df_nomatch["reads"], marker = '.', s = 2, color = "blue")
plt.xlabel("10X barcodes")
plt.ylabel("Anchor reads")
plt.title("Full (Barcodes)")
plt.minorticks_on()
plt.savefig("05_SABE_1172_UNHESMSV_barcodes_vs_anchreads_unmatched.png")
plt.close()

df_nomatch_cut = df_nomatch.loc[df_nomatch['barcode'] <= 2000]
df_nomatch_cut = df_nomatch_cut.loc[df_nomatch_cut['reads'] <= 100]

plt.scatter(df_nomatch_cut["barcode"], df_nomatch_cut["reads"], marker = '.', s = 2, color = "blue")
plt.xlabel("10X barcodes")
plt.ylabel("Anchor reads")
plt.title("Cutoff (Barcodes)")
plt.minorticks_on()
plt.savefig("05_SABE_1172_UNHESMSV_barcodes_vs_anchreads_cutoff.png")
plt.close()

plt.hist(df_nomatch_cut["barcode"], bins = 1000)
plt.title("Histogram of cutoff (Barcodes)")
plt.savefig("05_SABE_1172_UNHESMSV_barcodes_histogram_cutoff.png")
plt.close()

## 10X reads
plt.scatter(df_nomatch["x10reads"], df_nomatch["reads"], marker = '.', color = "blue")
plt.xlabel("10X reads")
plt.ylabel("Anchor reads")
plt.title("Full (10X reads)")
plt.minorticks_on()
plt.savefig("05_SABE_1172_UNHESMSV_10xreads_vs_anchreads_unmatched.png")
plt.close()

df_nomatch_10xreads_cut = df_nomatch.loc[df_nomatch['x10reads'] <= 100000]
df_nomatch_10xreads_cut = df_nomatch_10xreads_cut.loc[df_nomatch_10xreads_cut['reads'] <= 100]

plt.scatter(df_nomatch_10xreads_cut["x10reads"], df_nomatch_10xreads_cut["reads"], marker = '.', s = 2, color = "blue")
plt.xlabel("10X reads")
plt.ylabel("Anchor reads")
plt.title("Cutoff (10X reads)")
plt.minorticks_on()
plt.savefig("05_SABE_1172_UNHESMSV_10xreads_vs_anchreads_cutoff.png")
plt.close()


plt.hist(df_nomatch_10xreads_cut["x10reads"], bins = 1000)
plt.title("Histogram of cutoff (10X reads)")
plt.savefig("05_SABE_1172_UNHESMSV_10xreads_histogram_cutoff.png")
plt.close()

## with PMC
## Barcodes
df_nomatch_pmc = df_nomatch.loc[df_nomatch['nrs'].isin(supported_pmc)]

plt.scatter(df_nomatch["barcode"], df_nomatch["reads"], marker = '.', s = 2, color = "blue")
plt.scatter(df_nomatch_pmc["barcode"], df_nomatch_pmc["reads"], marker = 'x', s = 2, color = "purple")
plt.xlabel("10X barcodes")
plt.ylabel("Anchor reads")
plt.title("Full (Barcodes w/ PMC)")
plt.minorticks_on()
plt.savefig("05_SABE_1172_UNHESMSV_barcodes_vs_anchreads_unmatched_PMC.png")
plt.close()

df_nomatch_pmc_cut = df_nomatch_pmc.loc[df_nomatch_pmc['barcode'] <= 2000]
df_nomatch_pmc_cut = df_nomatch_pmc_cut.loc[df_nomatch_pmc_cut['reads'] <= 100]

plt.scatter(df_nomatch_cut["barcode"], df_nomatch_cut["reads"], marker = '.', s = 2, color = "blue")
plt.scatter(df_nomatch_pmc_cut["barcode"], df_nomatch_pmc_cut["reads"], marker = 'x', s = 2, color = "purple")
plt.xlabel("10X barcodes")
plt.ylabel("Anchor reads")
plt.title("Cuttoff (Barcodes w/ PMC)")
plt.minorticks_on()
plt.savefig("05_SABE_1172_UNHESMSV_barcodes_vs_anchreads_unmatched_PMC_cutoff.png")
plt.close()

## Barcodes histogram
import numpy as np

#Extract barcode values
bcs = list()
for ctg, val in X10_dict.items():
    if int(X10_dict[ctg]["bc"]) <= 5000:
        bcs.append(int(X10_dict[ctg]["bc"]))
bcs = np.array(bcs)

mu = np.mean(bcs)
sigma = np.std(bcs)
count, bins, ignored = plt.hist(bcs, 30, density = True)

plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (bins - mu)**2 / (2 * sigma**2) ), linewidth = 2, color = 'r')
plt.savefig("05_SABE_1172_UNHESMSV_barcodes_histogram.png")
plt.close()

####
print("\nAll done, have a nice day!")
