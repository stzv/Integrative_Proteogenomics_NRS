import collections
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from glob import glob
from math import log2
import math

import upsetplot


##
arms_info = "ARMS_bam_cluster_ProcessingResults/ARMS_Information.txt"
presto_info = "PRESTO_ProcessingResults/Presto_Information.txt"
thorax_info = "SABE_1172_UNHESMSV_RNASeq/THORAX_data/transcript_assembly/Thorax_Sample_Info.txt"

info_files = [arms_info, presto_info, thorax_info]
####
def Load_Info(tissue):
    SAMPLE_info = collections.defaultdict(dict)
    if tissue == "Presto":
        for line in open(presto_info, "r"):
            if line.startswith("#"): continue
            linesplit = line.split(";")
            UsedInRNASeq = linesplit[5]
            PatID = linesplit[12]
            Sex = linesplit[22]
            Age = linesplit[24]
            Smoker = linesplit[29]
            PackYears = linesplit[30]
            COPD = linesplit[32]
            if COPD == "0":
                Classification = "Control"
            else:                
                Classification = "COPD"
            ##
            if UsedInRNASeq == "Y":
                SAMPLE_info[PatID] = {tissue: PatID, "Info": (Sex, Age, Classification, Smoker)}
    elif tissue == "Thorax":
        for line in open(thorax_info, "r"):
            if line.startswith("#"): continue
            Key, PatID, COPD, Sex, Age, PackYears, FEV = line.strip("\n").split("\t")
            PatID =  "OE1507_" + PatID
            ## Smoker status
            if PackYears == "NA":
                continue
            else:
                Smoker = "0"
            ## Classification
            if COPD == "yes":
                Classification = "COPD"
            else:
                Classification = "Control"
            SAMPLE_info[PatID] = {tissue: PatID, "Info": (Sex, Age, Classification, Smoker)}
    else:
        for line in open(arms_info, "r"): 
            if line.startswith("#"): continue
            PatID, BiopID, NoseID, LungID, Sex, Age, Classification, Smoker, PackYears = line.strip("\n").split("\t")
            if Classification == "N/A": # Bad samples
                continue
            ## Smoker status
            if Smoker == "Never smoker":
                Smoker = "2"
            elif Smoker == "Past smoker":
                Smoker = "0"
            else:
                Smoker = "1"
            ## Sample IDs         
            if not BiopID == "N/A":
                BiopID = "104508-001-" + BiopID.rjust(3, '0')
            if not NoseID == "N/A":
                NoseID = "104508-001-" + NoseID.rjust(3, '0')
            if not LungID == "N/A":
                LungID = "104508-001-" + LungID.rjust(3, '0')
            ## Add to dict
            SAMPLE_info[PatID] = {"biopsy": BiopID, 
                                "nasalbrush": NoseID, 
                                "lungbrush": LungID, 
                                "Info": (Sex, Age, Classification, Smoker)}
    ##
    return SAMPLE_info

# Print simplified patient information table
outfile = open("Patient_Information_All_Cohorts.txt", "w")
outfile.write("#For ARMS samples, sampleIDs in order biopsy - nasalbrush - lungbrush\n")
outfile.write("# Smoking status: 0 ex smoker, 1 current smoker, 2 never smoker\n")
outfile.write("#PatID\tSex\tAge\tClassification\tSmoker\tCohort\tTissueID\n")
for infil in [arms_info, presto_info, thorax_info]:
    cohort = infil.split("/")[-1].split("_")[0]
    ## Information about samples
    SAMPLE_info = Load_Info(cohort)
    ##
    for patID, information in SAMPLE_info.items():
        Sex, Age, Classification, Smoker = information["Info"]
        if "ARMS" in patID:
            SampleID = "\t".join([information["biopsy"], information["nasalbrush"], information["lungbrush"]])
        else:
            SampleID = information[cohort]
        ##
        if any("#N/A" in i for i in (SampleID, Sex, Age, Classification)):
            continue
        outline = "\t".join([patID, Sex, Age, Classification, Smoker, cohort, SampleID])
        outfile.write(outline + "\n")

#
##
nrs_expression_arms = glob("ARMS_bam_cluster_ProcessingResults/03_NRS_coverage_freq_ARMS_*.txt")
nrs_expression_thorax = "SABE_1172_UNHESMSV_RNASeq/THORAX_data/04_SABE_1172_UNHESMSV_RNASeq_freq_Thorax.txt"
nrs_expression_presto = "PRESTO_ProcessingResults/03_NRS_coverage_freq_Presto.txt"

expression_files = nrs_expression_arms + [nrs_expression_presto]

def Get_Patient_Info(fil):
    NRS_PAT_info = collections.defaultdict(set)
    Cohort_Population = 0
    for line in open(fil, "r"):
        if line.startswith("#"): 
            Cohort_Population = len(line.strip("\n").split("\t")[1:])
            continue
        ##
        nrs, freq, samples = line.strip("\n").split("\t")
        ##
        for sample in samples.split(","):
            sample_id = sample.split(":")[0]
            PatID = [PatID for PatID, values in SAMPLE_info.items() if values[tissue] == sample_id]
            if PatID:
                info_sex = SAMPLE_info[PatID[0]]["Info"][0]
                info_disease = SAMPLE_info[PatID[0]]["Info"][2]
                NRS_PAT_info[nrs].add(",".join([PatID[0], info_sex, info_disease]))
    ##
    return NRS_PAT_info, Cohort_Population



df3 = pd.DataFrame()

for infil in  sorted(expression_files):
    ##
    NRS_phenotypes = collections.defaultdict(set)
    cohort_phenotypes = set()
    ##
    tissue = infil.split("_")[-1].replace(".txt", "")
    if tissue == "Presto":
        cohort = "Presto"
    else:
        cohort = " ".join(["ARMS", tissue.replace("brush", " brush")])
    #print(tissue, "\n")
    ## Information about samples
    SAMPLE_info = Load_Info(tissue)
    ## Add patient info to NRS, get cohort population size
    NRS_PAT_info, Cohort_Population = Get_Patient_Info(infil)
    ##
    for nrs, info in NRS_PAT_info.items():
        ##
        controls = 0
        for patient in info:
            NRS_phenotypes[nrs].add(patient)
            cohort_phenotypes.add(patient)
            if "Control" in patient:
                controls += 1
    ##
    all_pats = len(cohort_phenotypes) # All patients in the population
    all_control = len([c for c in cohort_phenotypes if "Control" in c]) # All control patients in population
    all_sick = len([c for c in cohort_phenotypes if not "Control" in c]) # All sick patient in population
    ##
    nrs_phenotypes = list()
    all_controls = len([c for c in cohort_phenotypes if "Control" in c])
    count = 0
    ##
    for nrs, phens in NRS_phenotypes.items():
        controls_WNRS = len([c for c in phens if "Control" in c]) # Controls with NRS
        PatNumWithNRS = len(phens) # All patients with NRS
        sick_WNRS = len([c for c in phens if not "Control" in c]) # Sick with NRS
        ##
        controls_WoNRS = all_control - controls_WNRS
        sick_WoNRS = all_sick - sick_WNRS
        ##
        odds_sick = (sick_WNRS + 0.5) / (sick_WoNRS + 0.5)
        odds_control = (controls_WNRS + 0.5) / (controls_WoNRS + 0.5)
        odds_ratio = odds_sick / odds_control
        ## Population frequency of NRS across all cohorts
        popul_freq = PatNumWithNRS / len(cohort_phenotypes)
        if not popul_freq >= 0.1: continue
        if log2(odds_ratio) > 0:
            count += 1
        # No need to filter popul freq, keeping all NRS which are >10% in own cohort
        nrs_phenotypes.append({"NRS": nrs,
                            "NRS_short": "_".join(nrs.split("_")[:2]),
                            "Population Frequency": round(popul_freq, 2),
                            "Sick_w_NRS": sick_WNRS,
                            "Control_w_NRS": controls_WNRS,
                            "Odds Ratio": round(odds_ratio, 2),
                            "Log2 Odds Ratio": round(log2(odds_ratio), 2),
                            "Samples": ";".join(phens),
                            "Cohort": cohort})

    df2 = pd.DataFrame(nrs_phenotypes).sort_values(by = "Population Frequency", ascending = False).reset_index(drop = True)
    df3 = df3.append(df2, ignore_index=True)
    ##
    # print("Total NRS >= 10% PopFreq", len(df2["NRS"]))
    # print("NRS count with OR > 0", count)
    # print("    Which is", round(count*100/len(nrs_phenotypes), 2), "%")
    # positive_odds = df2[df2["Log2 Odds Ratio"] > 0]
    # print("NRS with positive Odds Ratio", len(positive_odds["NRS"]))
    # # Correlation
    # correlation = df2["Population Frequency"].corr(df2["Log2 Odds Ratio"])
    # print("Correlation of Pop Freq with Log2 Odds:", round(correlation, 3))
    # print("\n")
    ## Save file
    df2.to_csv(f"01_NRS_ControlsProportion_{tissue}.tsv", sep = "\t", header = True, index = False)
    df3.to_csv(f"01_NRS_ControlsProportion.tsv", sep = "\t", header = True, index = False)
    ##



fig = plt.figure(figsize = (7, 7))
ax = fig.add_subplot()

ax = sns.histplot(data = df3, x = "Log2 Odds Ratio", hue = "Cohort",
             stat = "percent", bins = 30,
             palette = "tab10",
             alpha = 0.5,
             multiple = "stack", 
             element = "step",
             legend = False
             )

ax.axvline(0, linewidth = 0.5, dashes = [6, 2])
ax.set_xlabel("Log2 Odds Ratio")

plt.suptitle("Positive OR tendency across cohorts\nStacked histogram")

fig.tight_layout()
plt.savefig(f"01_NRS_OR.png")
