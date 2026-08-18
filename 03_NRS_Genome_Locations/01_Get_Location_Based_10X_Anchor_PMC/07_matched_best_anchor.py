import sys

print("Finding the best position for genome mapped NRS\n")

pmc_dict, anch_dict = dict(), dict()

## Create PMC library
for line in open("SABE_1172_UNHESMSV_genotyping/SABE_1172_UNHESMSV_PMC_locs.txt", "r"):
    ctg, pmc = line.strip("\n").split("\t")
    ctg = ctg.split("_")[0] + "_" + ctg.split("_")[1]
    pmc_dict[ctg] = pmc.split(":") #Split into flag, chr, pos, mapq, cigar

print("Total NRS with PMC", len(pmc_dict))

## Create anchor library
for line in open("05_SABE_1172_UNHESMSV_10x_anchor_match.txt", "r"):
    if not line.startswith("k141"):
        continue
    #
    ctg, pmc, x10, bc, reads, anchors_matched, anchors_all = line.strip("\n").split("\t")
    anchors_matched = anchors_matched.replace(f"{reads}|", "").split(",")
    anch_dict[ctg] = anchors_matched

print("Total NRS with Anchors", len(anch_dict), "\n")

## Choose the anchor pair for each NRS
for ctg in anch_dict:
    # Fetch PMC loc if exists
    #pmc = "N/A"
    #if ctg in pmc_dict.keys():
    #    pmc = pmc_dict[ctg]
    #
    if ctg in pmc_dict.keys():
        print(ctg, pmc_dict[ctg], anch_dict[ctg])
        sys.exit()
