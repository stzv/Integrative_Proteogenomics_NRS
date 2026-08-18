infile = open("01_SABE_UNHESMSV_1172_10x_merged_final2.txt", "r")

locs_reads_sum = dict()

print("Counting the reads and barcodes")

for line in infile:
    ctg = line.split("\t")[0]
    #print(ctg)
    hits = line.strip("\n").split("\t")[1].split(",")
    #
    locs = list()
    locs_sum = list()
    # Create list for the NRS
    if ctg not in locs_reads_sum:
        locs_reads_sum[ctg] = list()
    # Make list of chromosomes
    for h in hits:
        loc = h.split(":")[0]
        if loc not in locs:
            locs.append(loc)
    # Sum the reads for the locations
    clusters = list()
    for l in locs:
        subset = [h for h in hits if h.startswith(f"{l}:")]
        prev = None
        group = []
        for s in subset:
            if not prev or abs(int(s.split(":")[1]) - int(prev.split(":")[1])) <= 100:
                group.append(s)
            else:
                group = [s]
            prev = s
            if group not in clusters:
                clusters.append(group)
    # Process each cluster
    for cl in clusters:
        # Select top loc of cluster (highest # of own reads)
        cl_sorted = sorted(cl, key = lambda x: int(x.split(':')[2]), reverse = True)
        cl_top = ":".join(cl_sorted[0].split(":")[:3]) #keep chr, pos, own reads
        if ("_" in cl_top) and ("decoy" not in cl_top):
            #print(cl_top)
            cl_top = cl_top.replace(cl_top.split(":")[0], cl_top.split(":")[0].split("_")[1])
            #print(cl_top)
        # Count reads and barcodes for each cluster
        loc_reads = 0
        loc_bc = 0
        sample_count = list()
        for c in cl:
            loc_reads += int(c.split(":")[3])
            loc_bc += int(c.split(":")[4])
            sid = c.split(":")[5]
            if sid not in sample_count:
                sample_count.append(sid)
        if loc_bc > 1:
            locs_reads_sum[ctg].append(f"{cl_top}|{loc_reads}|{loc_bc}|{len(sample_count)}")
    #

#
outfile = open("02_SABE_UNHESMSV_1172_10Xlocs_readsbc_sum2.txt", "w+")

for ctg in locs_reads_sum:
    if locs_reads_sum[ctg]:
        # Sort locations from most barcodes to least
        locs_sorted = sorted(locs_reads_sum[ctg], key = lambda x: int(x.split('|')[2]), reverse = True)
        # Extract location with most barcodes as top position
        top_pos = locs_sorted[0]
        # Save into file
        outfile.write(f"{ctg}\t{top_pos}\t{','.join(locs_sorted)}\n")

outfile.close()


print("All done, have a nice day!")

