from glob import glob
import gzip

print("Load in NRS")
NRS_dict = dict()

for line in open("SABE_1172_UNHESMSV_genotyping/SABE1172_UNHESMSV_NRS_dark_freeze_final.fa", "r"):
    if line.startswith(">"):
        ctg = line.split("\t")[0].replace(">", "").split(" ")[0].strip("\n")
        length = int(line.split("_")[3].split(" ")[0]) - int(line.split("_")[2]) + 1
        NRS_dict[ctg] = {"Length": length}

print("Total", len(NRS_dict), "NRS")

print("Process 10X samples")
sample_files = glob("/mnt/fedot21/chromium/20210822_SABE_UNHESMSV_linkedreads_final2/*_locsSABE1172_UNHESMSV_NRS_dark_freeze_final2.txt.gz")
locs_dict = dict()
counter = 0

for sf in sample_files:
    sampleid = sf[sf.find("_linkedreads_final2/")+20:sf.find("_locsSABE1172")]
    #
    counter += 1
    print(counter, f"/{len(sample_files)}")
    # Open data
    data = gzip.open(sf, "r")
    next(data) # skip header line
    #
    locations = set()
    # Process data
    for line in data:
        Contig,Unique_barcodes,On_target_barcodes,On_target_reads,Total_reads,On_target_per,Hits = line.decode().replace("\n", "").split("\t")[:7]
        target = Hits.split('\t')[0]
        #print(Contig, f"{target}:{sampleid}:{On_target_barcodes}\n")
        if not Contig in locs_dict:
            locs_dict[Contig] = list()
        locs_dict[Contig].append(f"{target}:{On_target_reads}:{On_target_barcodes}:{sampleid}")
    #

locs = open("01_SABE_UNHESMSV_1172_10x_merged_final2.txt", "w+")

for n in locs_dict:
    if n in NRS_dict.keys():
        locs.write(f"{n}\t{','.join(locs_dict[n])}\n")
