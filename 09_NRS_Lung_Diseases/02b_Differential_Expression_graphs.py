import collections
from glob import glob
from math import log
import upsetplot
import matplotlib.pyplot as plt
import numpy as np
from venn import venn
import seaborn as sns

import sys

## Commented out just for the purpose of drawing the volcano plot at the end
# UP_dict = collections.defaultdict(dict)
# DOWN_dict = collections.defaultdict(dict)
# DE_dict = collections.defaultdict(dict)
# total = set()

# for infil in sorted(glob("02_DEGenes_*_COPDvsControl.txt")):
#     cohort = infil.replace("02_DEGenes_", "").replace("_COPDvsControl.txt", "").replace("_", " ")
#     cohort_list_up, cohort_list_down, cohort_list = list(), list(), list()
#     for line in open(infil, "r"):
#         # Skip header
#         if line.startswith("genes"): continue
#         #
#         gene, logFC, logCPM, F, PValue, FDR, Symbol, Description = line.strip("\n").split("\t")
#         if not gene.startswith("k141"): continue
#         cohort_list.append(gene)
#         total.add(gene)
#     DE_dict[cohort] = cohort_list
#     ##
#     print(cohort)
#     print(" Total DE:", len(cohort_list), ", from which NRS", len([x for x in cohort_list if x.startswith("k141")]))

# # total_nrs = [x for x in total if x.startswith("k141")]
# # print("\nTotal genes & NRS DE Expressed", len(total), ", from which NRS", len(total_nrs))

# ### Upset Plots
# def DrawUpsetPlot(dictionary, keyword):
#     data = upsetplot.from_contents(dictionary)
#     fig = plt.figure(figsize = (7, 6)) # width, height
#     ax = fig.add_subplot()
#     Upset_Intersections = upsetplot.UpSet(data,
#         sort_by = "cardinality", sort_categories_by = "input",
#         subset_size = "count", show_counts = True,
#         )
#     Upset_Intersections.plot()
#     plt.savefig(f"02b_UpsetPlot_{keyword}.png")
#     ##
#     return

# # DrawUpsetPlot(UP_dict, "UpRegulated")
# # DrawUpsetPlot(DOWN_dict, "DownRegulated")
# DrawUpsetPlot(DE_dict, "DE")


# ####
# healthy, sick = [], []
# for infil in glob("Patient_Information_ARMS*.txt") + ["Patient_Information_Presto.txt"]:
#     for line in open(infil, "r"):
#         if line.startswith("PatID"):
#             continue
#         patid = line.split("\t")[0]
#         classification = line.split("\t")[3]
#         if classification == "Control":
#             healthy.append(patid.replace("X", "").replace(".", "-"))
#         else:
#             sick.append(patid.replace("X", "").replace(".", "-"))

# #### Population frequency of DE NRS
# freq_input_files = glob("MOUNTED_FOLDERS/ARMS_bam_cluster_ProcessingResults/03_NRS_coverage_freq_ARMS_*.txt") + ["MOUNTED_FOLDERS/PRESTO_ProcessingResults/03_NRS_coverage_freq_Presto.txt"]

# FREQ_dict = dict()
# Expr_dict = collections.defaultdict(set)
# total_expressed_nrs = set()

# for infile in sorted(freq_input_files):
#     freqs = list()
#     cohort = infile.split("/")[-1].replace("03_NRS_coverage_freq_", "").replace(".txt", "").replace("_", " ")
#     if not cohort == "Presto":
#         coh = cohort.replace("brush", " brush")
#     else:
#         coh = "Presto"
#     FREQ_dict[cohort] = collections.defaultdict(dict)
#     ##
#     de_list = DE_dict[cohort]
#     ##
#     for line in open(infile, "r"):
#         if line.startswith("#"): 
#             sample_count = len(line.split("\t")[1:])
#             continue # skip header
#         nrs, freq, samples = line.split("\t")
#         nrs_short = "_".join(nrs.split("_")[:2])
#         Expr_dict[coh].add(nrs_short)
#         total_expressed_nrs.add(nrs_short)
#         if nrs in de_list:
#             FREQ_dict[cohort][nrs_short] = (str(round(float(freq), 2)), samples.split(","))

# cmap = sns.color_palette()

# fig = plt.figure(figsize = (7, 6)) # width, height
# ax = fig.add_subplot()
# venn(Expr_dict, cmap = cmap, legend_loc = "upper left")
# plt.savefig("02b_Venn_NRS_ExpressedInCohort.png")#, transparent=True)
# print("Total expressed NRS", len(total_expressed_nrs))


# ### Genome frequency in SABE
# counter_rare, counter_common = 0 , 0
# GEN_FREQ_dict = collections.defaultdict(str)
# for line in open("SABE_1172_UNHESMSV_genotyping/SABE_1172_UNHESMSV_frequency_nrs.txt", "r"):
#     if line.startswith("#"): continue
#     nrs, sample_count, freq = line.strip("\n").split("\t")
#     GEN_FREQ_dict[nrs] = freq
#     if float(freq) < 0.05:
#         counter_rare += 1
#     if float(freq) > 0.99:
#         counter_common += 1

# print("Rare NRS", counter_rare)
# print("Common NRS", counter_common)

# ### Genomic locations of DE NRS
# GENMAP_fil = "SABE_1172_UNHESMSV_genomemapping/02_NRS_GENMAP.tsv"
# GENMAP_dict = collections.defaultdict(dict)

# BP_dict = collections.defaultdict(dict)

# for line in open(GENMAP_fil, "r"):
#     if line.startswith("#"): continue
#     NRS, length, GENMAP, ConsensLoc, Method, Match, PMC, X10, Anchor, RNASeq = line.rstrip("\n").split("\t")
#     GENMAP_dict[NRS] = GENMAP
#     BP_dict[NRS] = length

# ##
# OUT_dict = dict()
# counter = 0

# for cohort, de_list in DE_dict.items():
#     for de in de_list:
#         # Get locations
#         if not de.startswith("k141"): continue # Keep only DE NRS
#         de = "_".join(de.split("_")[:2])
#         if de not in OUT_dict.keys():
#             OUT_dict[de] = {"Location": "", "Cohort": list()}
#         location = GENMAP_dict.get(de, "NA")
#         if location.startswith("chr"):
#             counter += 1
#         # Save to dict
#         OUT_dict[de]["Location"] = location.split(":")
#         OUT_dict[de]["Cohort"].append(cohort.replace("_", " "))
# print("Total mapped NRS", counter)


# ## GeneCards genes associiated with Asthma and COPD
# asthma = list()
# for line in open("asthma_genes.txt", "r"):
#     asthma.append(line.strip("\t"))

# copd = list()
# for line in open("COPD_genes.txt", "r"):
#     copd.append(line.strip("\t"))


# ## Mapping with chimeric reads
# CJ_mapping = collections.defaultdict()
# cj_infils = ["MOUNTED_FOLDERS/PRESTO_ProcessingResults/Chimeric_Reads/06_SABE_1172_UNHESMSV_NRS_two_best_CJs_Presto.txt",
#              "MOUNTED_FOLDERS/ARMS_bam_cluster_ProcessingResults/Chimeric_reads/06_SABE_1172_UNHESMSV_NRS_two_best_CJs_ARMS.txt"]

# for infil in cj_infils:
#     cohort = infil.split("_")[-1].replace(".txt", "")
#     CJ_mapping[cohort] = collections.defaultdict(str)
#     for line in open(infil, "r"):
#         ##NRS\tFirstLocation\tSecondLocation\tPopFreq_FirstLocation\tPopFreq_SecondLocation\tRNASeq_FirstLocation\tRNASeq_SecondLocation
#         NRS, CJ1, CJ2 = line.split("\t")[:3]
#         short_nrs = "_".join(NRS.split("_")[:2])
#         CJ_mapping[cohort][short_nrs] = f"{CJ1}|{CJ2}"

# ##
# GTF = dict()
# for line in open("GTF/Homo_sapiens.GRCh38.105_genesonly.gtf", "r"):
#     chrom, source, feature, start, end, score, strand, something, attributes = line.split("\t")
#     if chrom not in GTF.keys():
#         GTF[chrom] = list()
#     gene_id = [attr for attr in attributes.split(";") if "gene_id " in attr][0].replace("\"", "")
#     gene_sym = [attr for attr in attributes.split(";") if "gene_name" in attr]
#     if gene_sym:
#         gene_sym = gene_sym[0].replace("\"", "")
#     else:
#         gene_sym = "-"
#     ##
#     association = list()
#     if gene_sym in asthma:
#         association.append("Asthma")
#     if gene_sym in copd:
#         association.append("COPD")
#     GTF[chrom].append([gene_id.replace("gene_id ", ""), gene_sym.replace(" gene_name ", ""), chrom, start, end, ";".join(association)])

# ##
# def Check_GTF(location):
#     gene_loc = set()
#     if not location[0].startswith("chr"):
#         pass
#     else:
#         chrom, start, end, strand  = location
#         GTF_sub = GTF[chrom.replace("chr", "")]
#         for gene in GTF_sub:
#             gene_id, gene_sym, gene_chrom, gene_start, gene_end, association = gene
#             dist_start = (int(gene_start) - 1000) <= int(start)
#             dist_end = int(end) <= (int(gene_end) + 1000)
#             if dist_start and dist_end:
#                 gene_loc.add(f"{gene_id}:{gene_sym}:{association}")
#     ##
#     return gene_loc


# ##
# gwas_files = glob("GWAS_*_associations.csv")
# gwas_bed = open("GWAS_SNPs.bed", "w")

# for fil in gwas_files:
#     for line in open(fil, "r"):
#         # Skip header
#         if "Variant" in line:
#             continue
#         # Split line to variables
#         gwas_variant, gwas_gene, gwas_trait, gwas_study, gwas_loc = line.strip("\n").split(",")
#         # If variant was not mapped but it's location is in name
#         if "Mapping not available" in gwas_loc:
#             if gwas_variant.startswith("chr"):
#                 gwas_loc = gwas_variant.split("-")[0].replace("chr", "")
#             else:
#                 continue
#         ##
#         for l in gwas_loc.split("|"):
#             gwas_chr, gwas_pos = l.split(":")
#             gwas_bed.write("\t".join([gwas_chr.replace("chr", ""), gwas_pos, gwas_pos, gwas_variant.split("-")[0]]) + "\n")


# def Check_GWAS(location):
#     gwas = list()
#     if not location[0].startswith("chr"):
#         pass
#     else:
#         chrom, start, end, strand = location
#         for fil in gwas_files:
#             for line in open(fil, "r"):
#                 # Skip header
#                 if "Variant" in line:
#                     continue
#                 # Split line to variables
#                 gwas_variant, gwas_gene, gwas_trait, gwas_study, gwas_loc = line.strip("\n").split(",")
#                 # If variant was not mapped but it's location is in name
#                 if "Mapping not available" in gwas_loc:
#                     if gwas_variant.startswith("chr"):
#                         gwas_loc = gwas_variant.split("-")[0].replace("chr", "")
#                     else:
#                         continue
#                 # Split the gwas location
#                 for l in gwas_loc.split("|"):
#                     gwas_chr, gwas_pos = l.split(":")
#                     # Compare to the NRS location
#                     if chrom.replace("chr", "") == gwas_chr:
#                         if (abs(int(start) - int(gwas_pos)) <= 100000) or (abs(int(end) - int(gwas_pos)) <= 100000):
#                             distance = abs(int(start) - int(gwas_pos))
#                             gwas.append(":".join([gwas_gene, l, str(distance), gwas_variant, gwas_trait, gwas_study]))                
#     ##
#     return gwas

# ## Novel peptides
# NovelPeps_dict = collections.defaultdict(set)
# for line in open("MOUNTED_FOLDERS/XOmics_PublicData/04_NovelPeptides_Merged.tsv"):
#     if line.startswith("#"):
#         continue
#     #Peptide	NRS	GENMAP	GTF	DarkTrans	ProteinInfo	BestSpectrum_PerCohort	cohort
#     peptide, nrs = line.split("\t")[:2]
#     NovelPeps_dict[nrs].add(peptide)



# ##
# outfile = open("02b_DE_GENMAP.tsv", "w")
# outfile.write("#NRS\tLength\tGenome_Freq\tGENMAP\tGene\tCohorts\tRNASeq_Pop_Freq\tSick/TotalPat\tCJMap:CJ1|CJ2\tNovelPeptides\tGWAS_Hits:gSNP_ID-Distance\n")

# outfile2 = open("02b_DE_GENMAP.bed", "w")
# outfile3 = open("02b_GWAS_hits.bed", "w+")

# outfile4 = open("02b_DE_NRS_GWAS.tsv", "w")
# outfile4.write("#NRS\tGENMAP\tGTF\tDistance\tgSNP\tChr\tPos\tgGene\tgTrait\tgStudy\n")

# counter = 0
# gwas_out = set()

# for de in OUT_dict.keys():    
#     coh = OUT_dict[de]["Cohort"]
#     length = BP_dict[de]
#     gen_freq = GEN_FREQ_dict[de]
#     freq, sick_pats, cjmap = list(), list(), set()
#     for cohort in coh:
#         freq.append(FREQ_dict[cohort][de][0])
#         sick_pats.append(f'{str(len([pat for pat in FREQ_dict[cohort][de][1] if pat.split(":")[0] in sick]))}/{len(FREQ_dict[cohort][de][1])}')
#         cjmap.add(CJ_mapping[cohort.split()[0]].get(de, "NA"))
#     ## Fetch GENMAP and GTF info
#     loc = OUT_dict[de]["Location"]
#     gene = Check_GTF(loc)
#     gwas = Check_GWAS(loc)
#     gwas_ids = set()
#     if gwas:
#         for g in gwas:
#             gwas_gene, chrom, pos, distance, gwas_variant, gwas_trait, gwas_study = g.split(":")
#             # For the file with all DE NRS
#             gwas_ids.add(f"{gwas_variant}-{distance}")
#             # For the file with only DE NRS <-> GWAS hit
#             outfile4.write("\t".join([de, ":".join(loc),  ",".join(gene), distance, gwas_variant, chrom, pos, gwas_gene, gwas_trait, gwas_study]) + "\n")
#             # GWAS SNP bed
#             gwas_out.add(" ".join([chrom, pos, pos, gwas_variant.split("-")[0]]))
#         counter += 1
#     ## Novel peptides
#     peps = NovelPeps_dict.get(de, "NA")
#     ## Write out to file
#     # Tsv file with results
#     outline = "\t".join([de, 
#                          length, 
#                          gen_freq, 
#                          ":".join(loc), 
#                          ",".join(gene), 
#                          ",".join(coh), 
#                          ",".join(freq), 
#                          ",".join(sick_pats), 
#                          ",".join(cjmap), 
#                          peps, 
#                          ",".join(gwas_ids)])
#     outfile.write(outline + "\n")
#     # NRS bed
#     if loc[0].startswith("chr"):
#         outline2 = " ".join([loc[0].replace("chr", ""), loc[1], loc[2], de])
#         outfile2.write(outline2 + "\n")

# print("NRS w/ GWAS overlap (<=> 100kb)", counter)

# for gwas in gwas_out:
#     outfile3.write(gwas + "\n")

# # Venn of trait overlaps
# from venn import venn
# GWAS_overlap = collections.defaultdict(set)

# for fil in gwas_files:
#     f = fil.split("_")[1]
#     for line in open(fil, "r"):
#         gwas_variant = line.split(",")[0].split("-")[0]
#         GWAS_overlap[f].add(gwas_variant)


# fig = plt.figure(figsize = (7, 6)) # width, height
# ax = fig.add_subplot()

# venn(GWAS_overlap)

# plt.savefig("02b_Venn_GWAS_Trait_overlap.png", transparent=True)

##### Volcano plot for paper
import pandas as pd

infils = glob("02_Genes_ARMS_*_COPDvsControl.txt") + ["02_Genes_Presto_COPDvsControl.txt"]
df = pd.DataFrame()

for infil in sorted(infils):
    cohort = infil.replace("02_Genes_", "").replace("_COPDvsControl.txt", "")
    if not cohort == "Presto":
        cohort = cohort.replace("_", " ").replace("brush", " brush")
    dat = pd.read_csv(infil, sep = "\t")
    dat["Cohort"] = cohort
    dat = dat[dat["genes"].str.contains('k141')]
    ##
    df = pd.concat([df, dat], ignore_index=True)
    
down = df[(df['logFC'] < 0) & (df['FDR'] <= 0.05)]
up = df[(df['logFC'] > 0) & (df['FDR'] <= 0.05)]

print("Downregulated", len(down))
print("Upregulated", len(up))
print("Total expressed NRS", len(df))

df_sign = df[df['FDR'] <= 0.05]
df_sign.to_csv("02b_DE_NRS_significant_merged.txt", sep = "\t", index = False)

color_up = "#C47B4E"    
color_down = "#4C72B0"   
color_ns = "#858585"    

fig = plt.figure(figsize = (10, 6)) 

plt.scatter(x = df['logFC'], y = df['FDR'].apply(lambda x: -np.log10(x)),
            s = 10, 
            c = color_ns, alpha = 0.35,
            label="Not significant"
)

plt.scatter(x = up['logFC'],y = up['FDR'].apply(lambda x: -np.log10(x)),
            s = 14,
            c = color_up,
            label = "Up-regulated")

plt.scatter(x = down['logFC'], y = down['FDR'].apply(lambda x: -np.log10(x)), 
            s = 14,
            c = color_down,
            label = "Down-regulated")


plt.xlabel("log₂ fold change", fontsize=12, labelpad=8)
plt.ylabel("−log₁₀(FDR)", fontsize=12, labelpad=8)

plt.xlim(-8, 8)
plt.axhline(-np.log10(0.05), color = "grey", linestyle = "--", linewidth = 0.5)

plt.legend(frameon=True, fontsize=12, edgecolor = "#B0B0B0", markerscale=1.8)

#plt.suptitle("Differential expression of NRS in disease-specific datasets", fontsize=16)
plt.tight_layout()

plt.savefig("02b_DE_NRS.png")