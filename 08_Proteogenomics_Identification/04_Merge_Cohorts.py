import collections
from glob import glob

Novel_peptides = dict()
infiles = ["03_Peptides_NotFound_Ensembl_ARMS.txt", "03_Peptides_NotFound_Ensembl_Presto.txt"]

variables = ["NRS", "cohort", "MSSamples", "GENMAP", "GTF", "DarkTrans", "ClosestMatchProteinInfo", "BestHyperscore", "BestSpectrum_PerCohort", "HighestSpectrumFile", "SecondHyperscore", "SecondSpectrumID", "SecondSpectrumFile"]

def NestedDict():
    nested_dict = dict()
    for var in variables:
        nested_dict[var] = list()
    ##
    return nested_dict

for infil in infiles:
    cohort = infil.replace(".txt", "").split("_")[-1]
    for line in open(infil, "r"):
        if line.startswith("#"): 
            continue
        Peptide, NRS, GENMAP, GTF, DarkTrans, ClosestMatchProteinInfo, ProteinInfo, MatchInfo, BestHyperscore,  BestSpectrum_PerCohort, SecondHyperscore, SecondSpectrumID = line.strip("\n").split("\t")
        if Peptide not in Novel_peptides.keys():
            Novel_peptides[Peptide] = NestedDict()
        ##
        peplist = f"02_peptides_list_{cohort}.txt"
        MSSamples = ""
        for line in open(peplist, "r"):
            pep = line.split("\t")[0]
            if pep == Peptide:
                MSSamples = line.split("\t")[1]
        ##
        sample = BestSpectrum_PerCohort.split("_c")[0]
        HighestSpectrumFile = glob(f"Public_data/01_MSMS_Alignment_Mar*_2023/{sample}*{cohort}_psm.tsv")[0]
        ##
        if SecondSpectrumID == "NA":
            second_sample = "NA"
            SecondSpectrumFile = "NA"
        else:
            second_sample = SecondSpectrumID.split("_c")[0]
            SecondSpectrumFile = glob(f"Public_data/01_MSMS_Alignment_Mar*_2023/{second_sample}*{cohort}_psm.tsv")[0]
        ##
        for key in Novel_peptides[Peptide].keys(): 
                Novel_peptides[Peptide][key].append(locals()[key])


##
outfile = open("04_NovelPeptides_Merged.tsv", "w")
outfile.write("#Peptide\t" + "\t".join(variables) + "\n")

outfile2 = open("04_Novel_peptides_list_publicdata.txt", "w")

bed_out = open("04_NovelPeptides.bed", "w")
bed = set()
for peptide, results in Novel_peptides.items():
    outfile2.write(peptide + "\n")
    outline = list()
    if len(results["BestHyperscore"]) > 1:
        mscount = sum([int(samples) for samples in results["MSSamples"]])
        results["MSSamples"] = [str(mscount)]
        #
        transcripts = results["DarkTrans"]
        results["DarkTrans"] = [";".join(transcripts)]
        #
        cohorts = results["cohort"]
        results["cohort"] = [";".join(cohorts)]
        #
        hyperscores = [float(hs) for hs in results["BestHyperscore"]]
        idx = hyperscores.index(max(hyperscores))
        results["BestHyperscore"] = results["BestHyperscore"][idx]
        results["BestSpectrum_PerCohort"] = results["BestSpectrum_PerCohort"][idx]
        results["HighestSpectrumFile"] = results["HighestSpectrumFile"][idx]
    ##
    if len(results["SecondHyperscore"]) > 1:
        second_hyperscores = [float(hs) for hs in results["SecondHyperscore"] if not hs == "NA"]
        second_idx = second_hyperscores.index(max(second_hyperscores))
        results["SecondHyperscore"] = results["SecondHyperscore"][second_idx]
        results["SecondSpectrumID"] = results["SecondSpectrumID"][second_idx]
        results["SecondSpectrumFile"] = results["SecondSpectrumFile"][second_idx]
        ##
    for var in variables:
        r = results[var]
        outline.append("".join(r))
    ##
    outfile.write(f"{peptide}\t" + "\t".join(outline) + "\n") # PRinting to file
    ##
    loc = "".join(results["GENMAP"])
    nrs = "".join(results["NRS"])
    if loc.startswith("chr"):
        c = loc.split(":")[0].replace("chr", "")
        s = loc.split(":")[1]
        e = loc.split(":")[2]
        strand = loc.split(":")[3]
        bed.add(" ".join([c, s, e, nrs, strand]))

for b in bed:
    bed_out.write(b + "\n")