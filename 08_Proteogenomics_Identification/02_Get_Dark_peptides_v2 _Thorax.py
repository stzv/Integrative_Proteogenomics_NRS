#/usr/bin/python3

from glob import glob
import collections
import pastaq
from glob import glob
import os



#################
cohort = "Thorax"


#################
peptides_out_file = f"02_peptides_list_{cohort}.txt"
spectrum_files = f"Public_data/01_MSMS_Alignment_Mar27_2023/*{cohort}_psm.tsv" 

all_spectrum_files = glob(spectrum_files) 

header = list()
accession_dict = collections.defaultdict(dict)


print("Total samples to process", len(all_spectrum_files))
count = 0
peptides = collections.defaultdict(list)
spectrum_dict = collections.defaultdict(set)

for spectrum_file in all_spectrum_files:
    sampleid = spectrum_file.split("/")[-1].replace("_SwissProt_*_normalSearchMSFragger_*_psm.tsv", "")
    for line in open(spectrum_file, "r"):
        # Skip header
        if line.startswith("Spectrum"):
            continue
        #
        spectrum_id, fil_loc, peptide = line.split()[:3]
        score = line.split("\t")[17]
        nrs = line.split("\t")[28]
        #
        peptides[peptide].append(",".join([spectrum_id, score, nrs]))


pepout = open(peptides_out_file, "w+")
pepout.write("#Peptide\tMSSamples\tHighestHyperscore\tSpectrumID\tNRS\tAllMatches\n")
counter = 0

for peptide, info in peptides.items():
    # Sort by hyperscore
    sorted_info = sorted(peptides[peptide], key = lambda item: float(item.split(",")[1]), reverse = True)
    HighestHyperscore = sorted_info[0].split(",")[1]
    BestSpectrumID = sorted_info[0].split(",")[0]
    BestMatchGene = sorted_info[0].split(",")[2]
    #
    if BestMatchGene.startswith("k141") or BestMatchGene.startswith("TRX"):
        outline = [peptide, str(len(sorted_info)), HighestHyperscore, BestSpectrumID, BestMatchGene, "|".join(sorted_info)]
        pepout.write("\t".join(outline) + "\n")
        counter += 1

print("Total peptides", counter)