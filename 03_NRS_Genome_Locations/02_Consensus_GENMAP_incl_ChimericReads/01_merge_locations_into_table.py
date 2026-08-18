import sys

#### LOCATION FILES
NRS = "SABE_1172_UNHESMSV_genotyping/SABE_1172_UNHESMSV_NRS_dark_freeze_final.fa"
x10 = "SABE_1172_UNHESMSV_genomemapping/Based_on_reads_newest/02_SABE_UNHESMSV_1172_10Xlocs_readsbc_sum2.txt" 
anchors = "SABE_1172_UNHESMSV_genomemapping/Based_on_reads_newest/03_SABE_UNHESMSV_1172_anchors_reads_count.txt"
PMC = "SABE_1172_UNHESMSV_genotyping/SABE_1172_UNHESMSV_PMC_locs.txt"
#
RNA = "SABE_1172_UNHESMSV_RNASeq/05_RNASeq_chimeric_junctions_merged.txt" 
chm13 = "/dark/20210721_compare_asm/SABE_vs_CHM13_best_hits.txt" 

#### LOCATION DICTIONARIES
print("Load NRS")
NRS_list = dict()
for line in open(NRS, "r"):
    # Skip non-ID lines
    if not line.startswith(">k141"):
        continue
    #
    nrs_id = line.split(" ")[0].split("_")[0].replace(">", "") + "_" + line.split(" ")[0].split("_")[1]
    length = int(line.split(" ")[0].split("_")[3]) - int(line.split(" ")[0].split("_")[2]) + 1
    NRS_list[nrs_id] = length

print(f"  Total {len(NRS_list):,}")
##
print("Load 10X")
X10_dict = dict()
for line in open(x10, "r"):
    nrs_id, toploc, locs = line.strip("\n").split("\t")
    nrs_id = nrs_id.split("_")[0] + "_" + nrs_id.split("_")[1]
    X10_dict[nrs_id] = toploc

print(f"  Total {len(X10_dict):,}")
##
print("Load anchors")
ANCH_dict = dict()
for line in open(anchors, "r"):
    ctg, anchors = line.strip("\n").split("\t")
    anchors = anchors.strip("[]").replace("'", "").split(", ")
    ANCH_dict[ctg] = anchors

print(f"  Total {len(ANCH_dict):,}")
##
print("Load PMC")
PMC_dict = dict()
for line in open("SABE_1172_UNHESMSV_genotyping/SABE_1172_UNHESMSV_PMC_locs.txt", "r"):
    ctg, position = line.strip("\n").split("\t")
    ctg = ctg.split("_")[0] + "_" + ctg.split("_")[1]
    PMC_dict[ctg] = position

print(f"  Total {len(PMC_dict):,}")
##
# print("Load CJ anchoring")
# CJ_dict = dict()
# for line in open(RNA, "r"):
#     nrs, cjloc, read_id = line.strip("\n").split("\t")
#     #Save to dictionary
#     if nrs not in CJ_dict:
#         CJ_dict[nrs] = []
#     CJ_dict[nrs].append(cjloc)

# print(f"  Total {len(CJ_dict):,}")
# ##
# print("Load Chm13")
# chm13_dict = dict()
# for line in open(chm13, "r"):
#     linesplit = line.strip("\n").split("\t")
#     nrs = linesplit[0].split("_")[0] + "_" + linesplit[0].split("_")[1]
#     chrom = linesplit[3].replace("CHM13_", "")
#     #
#     if nrs not in chm13_dict.keys():
#         chm13_dict[nrs] = []
#     #
#     chm13_dict[nrs] = chrom
    
#print(f"  Total {len(chm13_dict):,}")

#### Combine into table
outfile = open("01_SABE_1172_UNHESMSV_locations_table.tsv", "w+")
outfile.write(f"NRS\tlen_bp\tPMC\t10X\tAnchors\n")

master_table = dict()
for nrs in NRS_list.keys():
    ## NRS length
    lens = NRS_list[nrs]
    ## Locations
    pmc_loc, x10_loc, anch_loc= "NA","NA","NA"
    # PMC
    if nrs in PMC_dict.keys():
        pmc_loc = PMC_dict[nrs]
    # 10X
    if nrs in X10_dict.keys():
        x10_loc = X10_dict[nrs]
    # Anchors
    if nrs in ANCH_dict.keys():
        anch_loc = ",".join(ANCH_dict[nrs])
    # RNA
    # if nrs in CJ_dict.keys():
    #     rna_loc = ",".join(CJ_dict[nrs])
    # # Chm13
    # if nrs in chm13_dict.keys():
    #     chm13_loc = chm13_dict[nrs]
    # Save
    outfile.write(f"{nrs}\t{lens}\t{pmc_loc}\t{x10_loc}\t{anch_loc}\n")#\t{rna_loc}\t{chm13_loc}\n")

# print("All done, have a nice day")
# import subprocess
# subprocess.call("email_stepanka.pl 01_merge_locations GENMAP", shell = True)