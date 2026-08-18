import itertools
from collections import Counter
import gzip

####
print("Load in anchors")
ANCHORS_dict = dict()

for line in gzip.open("chromium/20210408_get_anchors/04_SABE_1172_UNHESMSV_anchors.txt.gz", "rb"):
    ctg = line.decode().split("\t")[0]
    anch = line.decode().replace("\n", "").split("\t")[1:]
    ANCHORS_dict[ctg] = anch

####
print("Counting reads")
locs_reads_count = dict()

counter = 0

for n in ANCHORS_dict:
    counter += 1
    ## Start NRS analysis with clean slate
    # 
    anchor_list = ANCHORS_dict[n]
    chromosome_list = set()
    #
    if n not in locs_reads_count:
        locs_reads_count[n] = list()
    #
    for al in anchor_list: 
        k = al.split(":")[0]
        chromosome_list.add(k)
    #
    clusters = list()
    for cl in chromosome_list:
        # Extract all locs for this chromosome
        subset = [i for i in anchor_list if i.startswith(f"{cl}:")]
        # Separate locations into <>500 clusters 
        prev = None
        group = []
        for s in subset:
            if not prev or abs(int(s.split(":")[1]) - int(prev.split(":")[1])) <= 500:
                group.append(s)
            else:
                group = [s]
            prev = s
            if group not in clusters:
                clusters.append(group)
    # Process each cluster
    for cl in clusters:
        # Save the loc's values if reads => 5
        if len(cl) > 4:
            locs_reads_count[n].append(f"{len(cl)}|{'|'.join(cl)}")

###
print("Saving output")
outfile = open("03_SABE_UNHESMSV_1172_anchors_reads_count.txt", "w+")

for n in locs_reads_count:
    # save only NRS that have any anchor (>=5) at all
    if not len(locs_reads_count[n]) > 0:
        continue
    # Keep only top hits anchors
    locs_reads_count_sorted = sorted(locs_reads_count[n], key=lambda x: int(x.split('|')[0]), reverse = True) # Sort anchors by amount of reads
    anchor_scores = list()
    for anch in locs_reads_count_sorted:
        a = anch.split("|")[0] # Extract locs read count
        anchor_scores.append(a) # Append to list
    max_anchor_score = max(anchor_scores) # Find the highest read count
    anchs = list()
    for anch in locs_reads_count_sorted:
        if anch.startswith(f"{max_anchor_score}|"): # Only keep anchor locations with the highest number of reads; might be multiple
            anchs.append(anch)
    # Save into file
    outfile.write(f"{n}\t{anchs}\n")

outfile.close()
   
####
print("\nAll done. Have a good day!")
