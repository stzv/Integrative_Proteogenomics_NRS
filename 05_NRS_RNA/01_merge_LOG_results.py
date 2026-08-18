from glob import glob
import subprocess
import os

print("Fetch BAM files")

bamfiles = glob("brazil/20210817_STAR/results/*/*Log.final.out")
print(f" Total {len(bamfiles)} finished BAM files")


print("Check if BAM files are not corrupted")
counter = 0
infiles = list()

for bam in bamfiles:
    counter += 1
    if os.path.getsize(bam) > 0:
        error = subprocess.call(f"samtools quickcheck {bam.replace('Log.final.out', '.bam')}", shell = True)
        if not error:
            infiles.append(bam)
    # Counter
    if counter % 100 == 0:
        print(f" Done", counter)

print(" Passed", len(infiles))

# LOG file
print("Merge LOG files")
outfile = open("01_RNASeq_LOG_merged.txt", "w+")
outfile.write("#RNA_ID\tinput_reads_#\tavg_input_len\tuniq_map_#\tuniq_map_%\tuniq_map_avglen\tuniq_splices_annot_#\tmulti_loci_perc\tmulti_toomanyloci_perc\tunmap_toomanymismatch_perc\tunmap_tooshort_perc\tunmap_other_perc\tchimeric_reads_#\tchimeric_reads_perc\n")

# Save used infile files in order to have the same used in downstream pipeline
rna_infiles = open("01_SABE_UNHESMSV_1172_RNASeq_infile_list.txt", "w+")
for infil in infiles:
    rna_infiles.write(f"{infil}\n")
rna_infiles.close()

counter = 0

for infile in infiles:
    counter += 1
    #print(counter, "/", len(infiles))
    ## Values for extraction
    rna_id = infile.split("/")[-1].strip("Log.final.out")
    line_counter = 0
    #
    no_input,avg_input = "", ""
    # Unique reads
    uniq_map_no, uniq_map_perc, avg_map_len, no_splic_annot = "", "", "", ""
    # Multimapping reads
    multi_loci, multi_toomanyloci = "", ""
    # Unmapped reads
    toomany_mismatch, too_short, other, chimeric_reads_no, chimeric_reads_perc = "", "", "", "", ""
    ## Extract information
    for line in open(infile, "r"):
        ##
        line_counter += 1
        info = line.split(" |	")
        ##
        if line_counter == 6:
            no_input = info[1].strip("\n")
            continue
        elif line_counter == 7:
            avg_input = info[1].strip("\n")
            continue
        elif line_counter == 9:
            uniq_map_no = info[1].strip("\n")
            continue
        elif line_counter == 10:
            uniq_map_perc = info[1].strip("\n")
            continue
        elif line_counter == 11:
            avg_map_len = info[1].strip("\n")
            continue
        elif line_counter == 13:
            no_splic_annot = info[1].strip("\n")
            continue
        elif line_counter == 25:
            multi_loci = info[1].strip("\n")
            continue
        elif line_counter == 27:
            multi_toomanyloci = info[1].strip("\n")
            continue
        elif line_counter == 30:
            toomany_mismatch = info[1].strip("\n")
            continue
        elif line_counter == 32:
            too_short = info[1].strip("\n")
            continue
        elif line_counter == 34:
            other = info[1].strip("\n")
            continue
        elif line_counter == 36:
            chimeric_reads_no = info[1].strip("\n")
            continue
        elif line_counter == 37:
            chimeric_reads_perc = info[1].strip("\n")
            continue
    ##
    outfile.write(f"{rna_id}\t{no_input}\t{avg_input}\t{uniq_map_no}\t{uniq_map_perc}\t{avg_map_len}\t{no_splic_annot}\t{multi_loci}\t{multi_toomanyloci}\t{toomany_mismatch}\t{too_short}\t{other}\t{chimeric_reads_no}\t{chimeric_reads_perc}\n")

outfile.close()