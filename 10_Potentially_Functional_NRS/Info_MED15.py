
import collections
from glob import glob

nrs_list = set()
nrs_loc = collections.defaultdict(str)
nrs_bp = collections.defaultdict(str)
for line in open("NRS_inregion.txt", "r"):
    nrs = line.split()[0]
    loc = line.split()[2]
    bp = line.split()[1]
    ##
    nrs_list.add(nrs)
    nrs_loc[nrs] = loc
    nrs_bp[nrs] = bp


####
counter_rare, counter_common = 0 , 0
GEN_FREQ_dict = collections.defaultdict(str)
for line in open("SABE_1172_UNHESMSV_genotyping/SABE_1172_UNHESMSV_frequency_nrs.txt", "r"):
    if line.startswith("#"): continue
    nrs, sample_count, freq = line.strip("\n").split("\t")
    if nrs in nrs_list:
        GEN_FREQ_dict[nrs] = freq

####
freq_input_files = glob("ARMS_bam_cluster_ProcessingResults/03_NRS_coverage_freq_ARMS_*.txt") + ["PRESTO_ProcessingResults/03_NRS_coverage_freq_Presto.txt"]

FREQ_dict = dict()
cohorts = collections.defaultdict(list)

for infile in sorted(freq_input_files):
    freqs = list()
    cohort = infile.split("/")[-1].replace("03_NRS_coverage_freq_", "").replace(".txt", "").replace("_", " ")
    FREQ_dict[cohort] = collections.defaultdict(dict)
    ##
    for line in open(infile, "r"):
        if line.startswith("#"): 
            sample_count = len(line.split("\t")[1:])
            continue # skip header
        ##
        nrs, freq, samples = line.split("\t")
        nrs_short = "_".join(nrs.split("_")[:2])
        ##
        if nrs_short in nrs_list:
            FREQ_dict[cohort][nrs_short] = (str(round(float(freq), 2)), samples.split(","))
            cohorts[nrs_short].append(cohort)


####
healthy, sick = [], []
cohort_phenotypes = set()
for infil in glob("../Patient_Information_ARMS*.txt") + ["../Patient_Information_Presto.txt"]:
    for line in open(infil, "r"):
        if line.startswith("PatID"):
            continue
        patid = line.split("\t")[0]
        cohort_phenotypes.add(patid)
        classification = line.split("\t")[3]
        if classification == "Control":
            healthy.append(patid.replace("X", "").replace(".", "-"))
        else:
            sick.append(patid.replace("X", "").replace(".", "-"))

all_pats = len(cohort_phenotypes) # All patients in the population
all_control = len(healthy) # All control patients in population
all_sick = len(sick) # All sick patient in population

####
DE_dict = collections.defaultdict(list)

for line in open('../02b_DE_NRS_significant_merged.txt', "r"):
    if not line.startswith("k141"): continue
    #genes	logFC	logCPM	F	PValue	FDR	Cohort
    nrs, logfc, logcpm, f, pval, fdr, coh = line.strip("\n").split("\t")
    if nrs in nrs_list:
        DE_dict[nrs].append(f"{coh}|{logfc}|{fdr}")

####
outfile = open("00_NRS_MED15_inf.tsv", "w")
outfile.write("NRS\tLength\tGenFreq\tGENMAP\tExprInCohorts\tFreqInCohorts\tSick/TotalPatients\tOddsRatio\tDGE\n")

for nrs in nrs_list:
    gen_freq = GEN_FREQ_dict[nrs]
    coh = cohorts[nrs]
    length = nrs_bp[nrs]
    loc = nrs_loc[nrs]
    de = DE_dict[nrs]
    #
    freq, sick_pats, cjmap = list(), list(), set()
    sick_WNRS, controls_WNRS = 0, 0
    for cohort in coh:
        freq.append(FREQ_dict[cohort][nrs][0])
        sicks = len([pat for pat in FREQ_dict[cohort][nrs][1] if pat.split(":")[0] in sick])
        total_pats= len(FREQ_dict[cohort][nrs][1])
        sick_pats.append(f'{str(sicks)}/{total_pats}')
        ##
        sick_WNRS += sicks
        controls_WNRS += total_pats-sick_WNRS
    #
    controls_WoNRS = all_control - controls_WNRS
    sick_WoNRS = all_sick - sick_WNRS
    ##
    odds_sick = (sick_WNRS + 0.5) / (sick_WoNRS + 0.5)
    odds_control = (controls_WNRS + 0.5) / (controls_WoNRS + 0.5)
    odds_ratio = round(odds_sick / odds_control, 2)
    # Tsv file with results
    outline = "\t".join([nrs, 
                         length, 
                         gen_freq, 
                         loc,
                         ",".join(coh), 
                         ",".join(freq), 
                         ",".join(sick_pats),
                         str(odds_ratio),
                         ",".join(de)
                         ])
    outfile.write(outline + "\n")