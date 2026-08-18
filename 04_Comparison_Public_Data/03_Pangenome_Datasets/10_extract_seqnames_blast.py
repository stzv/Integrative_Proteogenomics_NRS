print("Processing NCBI_NR results")
blastin = "SABE_1172_UNHESMSV_NCBINR/SABE_1172_UNHESMSV_NCBINR.m6"
list_output = open(f"SABE_1172_UNHESMSV_NCBINR_aligned_list.txt", "w")
f = open(blastin, "r")
entries = []

for line in f:
    #qaccver saccver pident length mismatch gapopen btop qstart qend sstart send evalue bitscore species
    ctg, subjname, IM, align_len, mismatches, gaps, btop, qstart, qend, subjstart, subjend, evalue, bitscore= line.split("\t")
    str_length = ctg.split("_")
    length = (int(str_length[3]) - int(str_length[2]) + 1)
    complet = length / float(align_len)
    if (float(IM) >= 99) and (complet >= 0.95):
        entries.append(f'{ctg},{length}')
        
for item in set(entries):
    list_output.write(f"{item}\n")
            
list_output.close()
