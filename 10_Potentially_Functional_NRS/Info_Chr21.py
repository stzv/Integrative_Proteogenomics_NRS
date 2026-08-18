
import collections
from glob import glob
from Bio import SeqIO


####
nrs_list = set()
nrs_loc = collections.defaultdict(str)
nrs_bp = collections.defaultdict(str)

outfile_bed = open("NRS_BED.bed", "w")

for line in open("SABE_1172_UNHESMSV_genomemapping/02_NRS_GENMAP.tsv", "r"):
    if not line.startswith("k141"):
        continue
    nrs = line.split("\t")[0]
    loc = line.split("\t")[2]
    bp = line.split("\t")[1]
    ##
    if nrs == "k141_24537181":
        loc = "chr21:8214500:8215000:*"
    if not loc.startswith("chr") or loc.startswith("chrEBV"):
        continue
    chrom, start, end, strand = loc.split(":")
    if chrom == "chr21" and int(start) >= 8200000 and int(end) <= 8254000:
        nrs_list.add(nrs)
        nrs_loc[nrs] = loc
        nrs_bp[nrs] = bp
        outfile_bed.write(line)

####
counter_rare, counter_common = 0 , 0
GEN_FREQ_dict = collections.defaultdict(str)
for line in open("SABE_1172_UNHESMSV_genotyping/SABE_1172_UNHESMSV_frequency_nrs.txt", "r"):
    if line.startswith("#"): continue
    nrs, sample_count, freq = line.strip("\n").split("\t")
    if nrs in nrs_list:
        GEN_FREQ_dict[nrs] = freq


####
FREQ_Public_dict = collections.defaultdict(dict)
for line in open("SABE_1172_UNHESMSV_RNASeq/04_SABE_1172_UNHESMSV_RNASeq_freq.txt", "r"):
    if line.startswith("#"):
        continue
    ##
    nrs, freq, samples = line.split("\t")
    sample_count = len([s for s in samples.split(",")])
    nrs_short = "_".join(nrs.split("_")[:2])
    FREQ_Public_dict[nrs_short] = sample_count

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
PEP_dict = collections.defaultdict(set)
for line in open("XOmics_PublicData/04_NovelPeptides_Merged.tsv", "r"):
    if line.startswith("#"):
        continue
    pep = line.split("\t")[0]
    nrs = line.split("\t")[1]
    PEP_dict[nrs].add(pep)

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
    nrs_short = "_".join(nrs.split("_")[:2])
    if nrs_short in nrs_list:
        DE_dict[nrs_short].append(f"{coh}|logFC:{logfc}|FDR:{fdr}")


####
outfile = open("00_NRS_Chr21_inf.tsv", "w")
outfile.write("NRS\tLength\tGenFreq\tGENMAP\tNovelPeptides\tPubDatSampleCount\tExprInCohorts\tFreqInCohorts\tSick/TotalPatients\tOddsRatio\tDGE\n")

outfile_bed = open("00_Chr21_NRS.bed", "w")
outfil_bed_expr = open("Expressed_NRS_BED.bed", "w")

nrs_keep = set()

for nrs in nrs_list:
    gen_freq = GEN_FREQ_dict.get(nrs, "0")
    coh = cohorts.get(nrs, "")
    length = nrs_bp.get(nrs, "")
    loc = nrs_loc.get(nrs, "")
    de = DE_dict.get(nrs, "")
    pep = PEP_dict.get(nrs, "")
    pd_freq = FREQ_Public_dict.get(nrs, "0")
    # NRS to keep for pangenome insertions
    if de:
        nrs_keep.add(nrs)
    elif len(coh) == 4:
        nrs_keep.add(nrs)
    elif pep:
        nrs_keep.add(nrs)
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
                        ",".join(pep),
                        str(pd_freq),
                        ",".join(coh), 
                        ",".join(freq), 
                        ",".join(sick_pats),
                        str(odds_ratio),
                        ",".join(de)
                        ])
    outfile.write(outline + "\n")
    ##
    if loc.startswith("chr") and not loc.startswith("chrEBV"):
        chrom, start, end, strand = loc.split(":")
        # BED file - all NRS
        outfile_bed.write("\t".join([chrom, start, end, nrs]) + "\n")
        # BED file - expressed NRS in public
        outfil_bed_expr.write("\t".join([chrom, start, end, nrs]) + "\n")


####
fasta_out = open("00_NRS_Chr21.fasta", "w")
handle = SeqIO.parse("SABE_1172_UNHESMSV_genotyping/SABE_1172_UNHESMSV_NRS_dark_freeze_final.fa", format = "fasta")
for entry in handle:
    nrs = "_".join(str(entry.id).split("_")[:2])
    if nrs in nrs_keep:
        fasta_out.write(f">{entry.id}\n{str(entry.seq)}\n")