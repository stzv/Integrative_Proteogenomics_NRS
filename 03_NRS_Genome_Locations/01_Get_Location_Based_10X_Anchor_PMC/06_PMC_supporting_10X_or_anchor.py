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
    barcodes = int(top10x.split("|")[2]) # Extract position's barcode sum
    xreads = int(top10x.split("|")[1])
    top10x = top10x.split("|")[0].replace("chr", "")
    X10_dict[ctg] = {"top10x": top10x, "bc": int(barcodes), "reads": int(xreads)}

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
for line in open("/home/stepankaz/SABE_1172_UNHESMSV_genotyping/SABE_1172_UNHESMSV_PMC_locs.txt", "r"):
    ctg, position = line.strip("\n").split("\t")
    ctg = ctg.split("_")[0] + "_" + ctg.split("_")[1]
    chrom, pos = position.split(":")[1:3]
    PMC_dict[ctg] = {"chrom": chrom, "pos": pos}
    

#
print("\nCombine locations")

counter_10x_only, counter_anchor_only, counter_noloc = 0, 0, 0
barcodes, pmc_10x, anchreads, pmc_anch = list(), dict(), list(), dict()

anchonly = open("06_SABE_1172_UNHESMSV_anchors_only.txt", "w+")
x10only = open("06_SABE_1172_UNHESMSV_10x_only.txt", "w+")
nolocation = open("06_SABE_1172_UNHESMSV_noloc.txt", "w+")

for ctg in NRS_dict:
    # If nrs has only 10x position
    if (ctg in X10_dict.keys()) and (not ctg in ANCHORS_dict.keys()):
        # Counter
        counter_10x_only += 1
        # For graph
        barcodes.append(X10_dict[ctg]["bc"])
        # If PMC position, check if 10X supported
        pmc_loc = ["NA"]
        if ctg in PMC_dict.keys():       
            pmc_loc = PMC_dict[ctg]
            x10 = X10_dict[ctg]["top10x"]
            x10_chr, x10_pos = x10.split(":")[:2]
            x10_pos = int(x10_pos) * 1000
            pmc_chr = pmc_loc["chrom"]
            pmc_pos = pmc_loc["pos"]
            # check if <> 100kb
            if (pmc_chr == x10_chr) and (abs(int(pmc_pos) - int(x10_pos)) <= 100000):
                pmc_10x[ctg] = {"pmc":":".join(PMC_dict[ctg]), "x10":x10}
                pmc_loc = f"{pmc_chr}:{pmc_pos}"
        # Save into file
        x10only.write(f"{ctg}\t{pmc_loc}\t{X10_dict[ctg]['top10x']}|{X10_dict[ctg]['bc']}|{X10_dict[ctg]['reads']}\n")
    # Else if nrs has only anchor
    elif (not ctg in X10_dict.keys()) and (ctg in ANCHORS_dict.keys()):
        # Counter
        counter_anchor_only += 1
        # For graph
        reads = list()
        for a in ANCHORS_dict[ctg]:
            reads.append(int(a.split("|")[0]))
        anchreads.append(max(reads))
        # If PMC position, check if any anchor <> 500bp
        pmc_loc = ["NA"]
        if ctg in PMC_dict.keys():
            pmc_chr = PMC_dict[ctg]["chrom"]
            pmc_pos = PMC_dict[ctg]["pos"]
            for anch in ANCHORS_dict[ctg]: # for each loc cluster
                match = "NO"
                for an in anch: # for each hit in a loc cluster
                    anc_chr = anch.split("|")[1].split(":")[0]
                    anc_pos = anch.split("|")[1].split(":")[1]
                    if (pmc_chr == anc_chr) and (abs(int(pmc_pos) - int(anc_pos)) <= 500):
                        match = "YES"
                if match == "YES":
                    pmc_anch[ctg] = {"pmc": ":".join(PMC_dict[ctg]),"anchor": anch}
                    pmc_loc = f"{pmc_chr}:{pmc_pos}"
        # Save into file
        anchonly.write(f"{ctg}\t{pmc_loc}\t{ANCHORS_dict[ctg]}\n")
    # If nrs does not have any location
    elif (not ctg in X10_dict.keys()) and (not ctg in ANCHORS_dict.keys()):
        counter_noloc += 1
        pmc_loc = ["NA"]
        if ctg in PMC_dict.keys():
            pmc_loc = PMC_dict[ctg]
        nolocation.write(f"{ctg}\t{pmc_loc}\n")
    #

#
print("NRS with only 10X", counter_10x_only)
print("Supported by PMC", len(pmc_10x))
print("NRS with only Anchor", counter_anchor_only)
print("Supported by PMC", len(pmc_anch))

# Plots
import pyplotlib as plt


## Barcodes
plt.scatter(range(0, len(barcodes)), barcodes, marker = '.', s = 2, color = "blue")
plt.ylabel("10X barcodes")
#plt.title("Full (Barcodes)")
plt.minorticks_on()
plt.savefig("06_SABE_1172_UNHESMSV_barcodes_10Xonly.png")
plt.close()

## Anchor reads
plt.scatter(range(0, len(anchreads)), anchreads, marker = '.', s = 2, color = "blue")
plt.ylabel("Anchor reads")
#plt.title("Full (Barcodes)")
plt.minorticks_on()
plt.savefig("06_SABE_1172_UNHESMSV_reads_Anchorsonly.png")
plt.close()
