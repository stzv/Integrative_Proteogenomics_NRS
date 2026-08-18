import subprocess
import os
from glob import glob
import numpy

import sys
from numpy.core.fromnumeric import transpose

#####
print("Load in NRS library")
nrs_library = dict()
NRS_fasta = "SABE_1172_UNHESMSV_genotyping/SABE_1172_UNHESMSV_NRS_dark_freeze_final.fa"

for line in open(NRS_fasta, "r"):
    if line.startswith(">"):
        ctg = line.strip("\n").split(" ")[0].strip(">")
        #length = int(line.strip("\n").split(" ")[0].split("_")[3]) - int(line.strip("\n").split("_")[2]) + 1
        #NRS_dict[ctg] = length
        nrs_library[ctg] = {"PopulFreq": 0, "TransFreq": 0, "Exon": set(), "RNAids": None, "Tissues": None, "CellTypes": None}


#####
print("Download SDRF files")

accessids_infile = open("03_SABE_1172_UNHESMSV_ACCESSIDs_RNAIDs.txt", "r")
accessids = dict()

for line in accessids_infile:
    domain = line.split(",")[0]
    accessids[domain] = line.strip("\n").split(",")[2:]
    if os.path.isfile(f"SABE_1172_UNHESMSV_RNASeq/SDRF_files/{domain}.sdrf.txt") == False:
        subprocess.call(f"wget -P SDRF_files https://www.ebi.ac.uk/arrayexpress/files/{domain}/{domain}.sdrf.txt", shell = True)


#####
print("Extract tissue information")

def GetInformation():
    accessinf_outfile = open(accessinf_file, "w+")
    accessinf_outfile.write("#Domain,RNASample,Organism,Tissue,CellType,Disease\n")
    ##
    sdrf_files = glob("SDRF_files/*.sdrf.txt")
    print(f" Total {len(sdrf_files)} SDRF files")
    ##
    sample_information, domain_dict = list(), dict()
    for fil in sdrf_files:
        # Line counter
        counter = 0
        # Environment preparation
        processed_domain = fil.strip(".sdrf.txt").split("/")[-1]
        processed_samples = accessids[processed_domain]
        sample_index, organism_index, tissue_index, cell_type_index, disease_index = [""] * 5
        # Get the information
        for line in open(fil, "r"):
            counter += 1
            line_split = line.split("\t")
            if counter == 1: # Header -> Find indexes of the information
                # Run name
                run_index = [ls for ls in line_split if "[ENA_RUN]" in ls]
                run_index = line_split.index(run_index[0])
                # Source name or ENA sample
                if any("[ENA_SAMPLE]" in ls for ls in line_split):
                    sample_index = [ls for ls in line_split if "[ENA_SAMPLE]" in ls]
                    sample_index = line_split.index(sample_index[0])
                elif any("Source Name" in ls for ls in line_split):
                    sample_index = [ls for ls in line_split if "Source Name" in ls]
                    sample_index = line_split.index(sample_index[0])
                else:
                    print("NO SOURCE NAME")
                    print(line_split)
                    sys.exit()
                # Organism
                organism_index = [ls for ls in line_split if ( ("Characteristics [organism]" in ls) or ("Characteristics[organism]" in ls) or ("Characteristics[Organism]" in ls) )]
                organism_index = line_split.index(organism_index[0])
                # Disease
                if any("[disease]" in ls for ls in line_split):
                    disease_index = [ls for ls in line_split if ("[disease]") in ls]
                    disease_index = line_split.index(disease_index[0])
                elif any("phenotype" in ls for ls in line_split):
                    disease_index = [ls for ls in line_split if ("phenotype") in ls]
                    disease_index = line_split.index(disease_index[0])
                elif any("infect" in ls for ls in line_split):
                    disease_index = [ls for ls in line_split if ("infect") in ls]
                    disease_index = line_split.index(disease_index[0])
                # Organism part or tissue
                if any("[tissue]" in ls for ls in line_split):
                    tissue_index = [ls for ls in line_split if ("[tissue]") in ls]
                    tissue_index = line_split.index(tissue_index[0])
                elif any("organism part" in ls for ls in line_split):
                    tissue_index = [ls for ls in line_split if ("organism part") in ls]
                    tissue_index = line_split.index(tissue_index[0])
                elif any("metastatic site" in ls for ls in line_split):
                    tissue_index = [ls for ls in line_split if ("metastatic site") in ls]
                    tissue_index = line_split.index(tissue_index[0])
                # Cell type
                if any("[cell type]" in ls for ls in line_split):
                    cell_type_index = [ls for ls in line_split if "[cell type]" in ls]
                    cell_type_index = line_split.index(cell_type_index[0])
                elif any("Sample_source_name" in ls for ls in line_split):
                    cell_type_index = [ls for ls in line_split if "[Sample_source_name]" in ls]
                    cell_type_index = line_split.index(cell_type_index[0])
                elif any("cell line" in ls for ls in line_split):
                    cell_type_index = [ls for ls in line_split if "[cell line]" in ls]
                    cell_type_index = line_split.index(cell_type_index[0])
                elif any("[strain]" in ls for ls in line_split):
                    cell_type_index = [ls for ls in line_split if "[strain]" in ls]
                    cell_type_index = line_split.index(cell_type_index[0])
                elif any("strain or line" in ls for ls in line_split):
                    cell_type_index = [ls for ls in line_split if "[strain or line]" in ls]
                    cell_type_index = line_split.index(cell_type_index[0])
                #strain or line
            else: # Process the samples
                # Only save samples that were actually processed
                sample = line_split[run_index]
                if not sample in processed_samples:
                    continue
                ## Use index for information retrieval
                if organism_index == "":
                    organism = "N/A"
                else:
                    organism = line_split[organism_index]
                #
                if tissue_index == "":
                    tissue = "N/A"
                else:
                    tissue = line_split[tissue_index]
                #
                if cell_type_index == "":
                    cell_type = "N/A"
                else:
                    cell_type = line_split[cell_type_index]
                #
                if disease_index:
                    disease = line_split[disease_index]
                else:
                    disease = "N/A"
                # Save the information
                result = f"{processed_domain};{sample};{organism};{tissue};{cell_type};{disease}\n"
                if result not in sample_information:
                    sample_information.append(result)
                    if sample not in domain_dict.keys():
                        domain_dict[sample] = {}
                    domain_dict[sample] = {"Organism": organism, "Tissue": tissue, "Cell type": cell_type, "Disease": disease}
    ##
    for i in sample_information:
        accessinf_outfile.write(i)
    ##
    return


accessinf_file = "07_SABE_1172_UNHESMSV_RNASeq_sample_informations_v2.txt"
if os.path.isfile(accessinf_file) == False:
    GetInformation()
else:
   domain_dict = dict()
   for line in open(accessinf_file, "r"):
      if line.startswith("#"): 
        continue
      ##
      Domain, RNASample, Organism, Tissue, CellType, Disease = line.strip("\n").split(";")
      # If tissue not applicable or not available - removed from the study, skip it for future steps
      if any(x == y for x in [Tissue, CellType] for y in ["not available", "not applicable"]):
        continue
      if RNASample not in domain_dict.keys():
         domain_dict[RNASample] = {}
      domain_dict[RNASample] = {"Organism": Organism, "Tissue": Tissue, "Cell type": CellType, "Disease": Disease}

#####
print("Load in NRS frequency in population")
pop_freq = {}

for line in open("SABE_1172_UNHESMSV_genotyping/SABE_1172_UNHESMSV_frequency_nrs.txt", "r"):
    if line.startswith("#"):
        continue
    nrs, sample_count, freq = line.strip("\n").split("\t")
    pop_freq[nrs] = freq


for nrs in nrs_library:
    nrs_short = "_".join(nrs.split("_")[:2])
    nrs_library[nrs]["PopulFreq"] = pop_freq[nrs_short]

del(pop_freq)

#####
print("Load in NRS transcription frequency")
rna_freq_file = "04_SABE_1172_UNHESMSV_RNASeq_freq.txt"

delim = "\t"

for line in open(rna_freq_file, "r"):
  line = line.strip("\n")
  if line.startswith("#"):
    continue
  #
  nrs_id = line.split("\t")[0]
  freq = line.split("\t")[1]
  rna_ids = line.split("\t")[2].split(",")
  nrs_library[nrs_id]["TransFreq"] = freq
  nrs_library[nrs_id]["RNAids"] = rna_ids


##### 
print("Load uniquely mapped reads information")
logout = "01_RNASeq_LOG_merged.txt"
rna_uqm = dict()

for line in open(logout, "r"):
  # Skip header
  if line.startswith("RNA_ID"):
    continue
  #
  rna_id = line.split("\t")[0]
  uniq_mapped = line.split("\t")[3]
  rna_uqm[rna_id] = uniq_mapped


##### 
print("Include Gene annotation")
geneannotation = open("06_SABE_1172_UNHESMSV_NRS_CJ_GTF.txt", "r")


for line in geneannotation:
  # Skip header
  if not line.startswith("k141"):
    continue
  #
  nrs, cj, gene_pos = line.strip("\n").split("\t")
  exon = gene_pos.split(":")[1]
  nrs_library[nrs]["Exon"].add(exon)

#####
print("Save it into final table")
final_out = open("07_SABE_1172_UNHESMSV_NRS_RNA_tissue_table.txt", "w+")
final_out.write("#NRS\tPopul_Freq\tRNA_Freq\tTissues\tExons\n")

cpm_out = open("07_SABE_1172_UNHESMSV_NRS_RNA_tissue_table_cpm>10.txt", "w+")

selection_out = open("07_SABE_1172_UNHESMSV_PublicRNA_selection.txt", "w+")
selection_out.write("#NRS;population_freq;transcription_freq;RNAsample;CPM;tissue\n")

newline = '\n'

for nrs, info in nrs_library.items():
  #
  population_freq = info["PopulFreq"]
  transcript_freq = info["TransFreq"]
  #
  expressed_gene = ",".join(info["Exon"])
  # Fetch tissues for each RNA sample the NRS is transcribed in
  RNAids_list = list()
  rnaids = info.get("RNAids", None)
  if not rnaids == None:
    rna_list = [i.split(":")[0] for i in info["RNAids"]]
    for rna, values in domain_dict.items():
      if rna in rna_list:
        #
        tissue = domain_dict[rna]["Tissue"].replace(newline, ' ').replace(",","_")
        if tissue == "N/A" or tissue == "n/a" or tissue == "":
          tissue = domain_dict[rna]["Cell type"]
        # Calculate CPM for each RNA Sample
        cpm = round(([int(i.split(":")[1]) for i in info["RNAids"] if i.split(":")[0] == rna][0] / int(rna_uqm[rna])) * 1000000, 2)
        #
        RNAids_list.append(f"{rna}:{cpm}:{tissue}")
    RNAids_list = sorted(RNAids_list, key=lambda x: x.rsplit(':')[1], reverse=True)
  # Print into file
  final_out.write(f"{nrs}\t{population_freq}\t{transcript_freq}\t{','.join(RNAids_list)}\t{expressed_gene}\n")
  # Print out well expressed contigs
  if (len(RNAids_list) > 0) and (int(float(RNAids_list[0].split(":")[1])) >= 10):
      cpm_out.write(f"{nrs}\t{population_freq}\t{transcript_freq}\t{','.join(RNAids_list)}\t{expressed_gene}\n")
  # Print out selection for validation - 50% population validation, > 0 RNA expression frequency, only highest expressed RNA:cpm:tissue
  if float(population_freq) >= 0.5 and (float(transcript_freq) > 0) and (len(RNAids_list) > 0):
    selected_tissue = next(iter(RNAids_list), "NO TISSUE INFORMATION")
    selection_out.write(f"{nrs};{population_freq};{transcript_freq};{selected_tissue};{expressed_gene}\n")


# ####
# print("Create graphs")

# import pandas as pd
# import matplotlib.pyplot as plt
# import matplotlib
# matplotlib.use("Agg")
# matplotlib.rcParams['agg.path.chunksize'] = 500000
# import numpy as np

# fig, ax = plt.subplots()

# data = transpose(pd.DataFrame(nrs_library))
# data["PopulFreq"] = pd.to_numeric(data["PopulFreq"])
# data["TransFreq"] = pd.to_numeric(data["TransFreq"])

# data.plot(x = "PopulFreq", y = "TransFreq", kind = "scatter")

# #ax.plot(pf, ef, ".")
# ax.set_xlabel("Population frequency")
# ax.set_ylabel("Transcription frequency")

# fig.tight_layout()
# plt.savefig("07_SABE_1172_UNHESMSV_populfreq_vs_transfreq.png")

# fig2, ax2 = plt.subplots()

# ax2.plot(pf, ef, ".")

# ax2.set_xlabel("Population frequency")
# ax2.set_ylabel("Transcription frequency")

# fig2.tight_layout()
# plt.savefig("07_SABE_1172_UNHESMSV_populfreq_vs_transfreq_cutoff.png")

# fig3, ax3 = plt.subplots()

# ax3.plot(pf, cpm_graph, ".")

# ax3.set_xlabel("Population frequency")
# ax3.set_ylabel("CPM")
# ax3.set_title("CPM of highest transcribed RNA sample")

# fig3.tight_layout()
# plt.savefig("07_SABE_1172_UNHESMSV_populfreq_vs_cpm_cutoff.png")
# ####
# subprocess.call("email_stepanka.pl NRS RNA table", shell = True)
# print("\nAll done, have a nice day!")