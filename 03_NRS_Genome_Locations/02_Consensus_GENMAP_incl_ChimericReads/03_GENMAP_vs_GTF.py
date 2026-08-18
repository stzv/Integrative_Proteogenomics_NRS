import collections

print("Load GTF")
GTF_dict = dict()

for line in open("GTF/Homo_sapiens.GRCh38.105_original.gtf", "r"):
    if line.startswith("#"): continue
    chrom, source, feature, start, end, score, strand, something, attributes = line.split("\t")
    if chrom not in GTF_dict.keys():
        GTF_dict[chrom] = list()
    if not feature == "gene":
        continue
    ##
    gene_id = [attr for attr in attributes.split(";") if "gene_id" in attr][0].replace("gene_id", "").strip(" \"")
    biotyp = [attr for attr in attributes.split(";") if "gene_biotype" in attr][0].replace("gene_biotype", "").strip(" \"")
    if "gene_name" in attributes:
        gene_sym = [attr for attr in attributes.split(";") if "gene_name" in attr][0].replace("gene_name", "").strip(" \"")
    else:
        gene_sym = ""
    GTF_dict[chrom].append([gene_id, gene_sym, chrom, start, end, feature, biotyp])


print("Process")
outfile = open("03_NRS_GENMAP_GTF.tsv", "w")
outfile_bed = open("03_NRS_GENMAP.bed", "w")
counter = 0
counter_gtf, counter_proteincoding = 0, 0

for line in open("02_NRS_GENMAP.tsv", "r"):
    if line.startswith("#"): continue
    NRS, Length_bp, genmap, Consens_Loc, GENMAP_Method, Matches, PMC, X10, Anchors, RNASeq = line.split("\t")
    counter += 1
    if not genmap.startswith("chr") or genmap.startswith("chrEBV"):
        GTF = "NA"
        outfile.write("\t".join([NRS, Length_bp, genmap, GTF, GENMAP_Method]) + "\n")
    else:        
        chrom, start, end, strand = genmap.replace("chr", "").split(":")
        gtf_subset = GTF_dict.get(chrom)
        GTF = set()
        for entry in gtf_subset:
            gene_id, gene_sym, gene_chrom, gene_start, gene_end, feature, biotyp = entry
            if int(gene_start) <= int(start) and int(gene_end) >= int(end): # Within gene
                GTF.add(":".join([gene_id, gene_sym, biotyp, gene_chrom, start, end]))
            elif abs(int(gene_start) - int(start)) <= 10000: # in front of gene
                GTF.add(":".join([gene_id, gene_sym, biotyp, gene_chrom, start, end]))
        if not GTF:
            GTF = "NA"
        else:
            GTF = ",".join(GTF)
        ##
        outfile.write("\t".join([NRS, Length_bp, genmap, GTF, GENMAP_Method]) + "\n")
        outfile_bed.write(" ".join([chrom, start, end, NRS]) + "\n")
    ##
    if not GTF == "NA":
        counter_gtf += 1
        if "protein_coding" in GTF:
            counter_proteincoding += 1
    if counter*100/600000 % 10 == 0:
        print(round(counter*100/600000, 2), "%")

print("Total NRS with GTF", counter_gtf)
print("In protein coding region", counter_proteincoding)