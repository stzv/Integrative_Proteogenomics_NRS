import os
import subprocess
from Bio import SeqIO
import gzip
import collections
import re
from pyopenms import ProteaseDB, ProteaseDigestion

###########
dataset = "Presto"

###########

peptides_list_file = f"02_peptides_list_{dataset}.txt"
novel_pep_outfile = f"03_Peptides_NotFound_Ensembl_{dataset}.txt"

####
print("Load in GENMAP")
GENMAP_fil = "03_NRS_GENMAP_GTF.tsv"
GENMAP_dict = collections.defaultdict(str)

for line in open(GENMAP_fil, "r"):
    if line.startswith("#"): continue
    NRS, Length_bp, genmap, GTF, GENMAP_Method = line.split("\t")
    GENMAP_dict[NRS] = "\t".join([genmap, GTF])

####
print("Load in GTF")

####
print("Load in Ensemble")
ensembl_fil = "ensembl_hs_pep/Homo_sapiens.GRCh38.pep.all.fa.gz"
swisprot_fil = "swisprot/Swisprot_true_format_fasta_includeIsoform_tr-2023.03.24-14.17.44.13.fasta"

ensembl = SeqIO.parse(gzip.open(ensembl_fil, "rt"),'fasta')
swisprot = SeqIO.parse(swisprot_fil,'fasta')

ensembl_seqs = collections.defaultdict(dict)
ensembl_ids = collections.defaultdict(str)

def GetEnsemblInfo(attributes):
    if "gene:" in attributes:
        gene_id = [attr for attr in attributes.split() if "gene:" in attr][0].replace("gene:", "")
    else:
        gene_id = ""
    if "gene_symbol" in attributes:
        gene_sym = [attr for attr in attributes.split() if "gene_symbol:" in attr][0].replace("gene_symbol:", "")
    else:
        gene_sym = ""
    if "gene_biotype:" in attributes:
        biotyp = [attr for attr in attributes.split() if "gene_biotype:" in attr][0].replace("gene_biotype:", "")
    else:
        biotyp = ""
    if ":GRCh38:" in attributes:
        location = [attr for attr in attributes.split() if ":GRCh38:" in attr][0]
    else:
        location = ""
    ##
    return gene_id, gene_sym, biotyp, location

for en in ensembl:
    ensembl_seqs[en.id] = str(en.seq)
    ##
    gene_id, gene_sym, biotyp, location = GetEnsemblInfo(str(en.description))
    ensembl_ids[en.id] = ",".join([gene_id, gene_sym, biotyp, location])

for en in swisprot:
    ensembl_seqs[en.id] = str(en.seq)
    ##
    gene_id, gene_sym, biotyp, location = GetEnsemblInfo(str(en.description))
    ensembl_ids[en.id] = ",".join([gene_id, gene_sym, biotyp, location])

####
print("Load in peptides information")

PEP_dict = collections.defaultdict(dict)
pep_set = set()
header = ""
for line in open(peptides_list_file, "r"):
    if line.startswith("#"):
        header = line
        continue
    #Peptide	MSSamples	HighestHyperscore	SpectrumID  SecondHyperscore    SecondSpectrumID	NRS	AllMatches
    Peptide = line.split("\t")[0]
    PEP_dict[Peptide] = line.split("\t")
    pep_set.add(Peptide)

####
if "Thorax" in dataset:
    TRXIDs = collections.defaultdict(str)
    for line in open("Thorax/08a_Thorax_SABE1172_UNHESMSV_IDs.txt"):
        trxid, nrsid = line.strip("\n").split("\t")
        TRXIDs[trxid] = nrsid

####
print(" Check for known peptides")
print("  Count to check", len(pep_set))

user_enzyme = "Trypsin/P"
insilico_digestion_object = ProteaseDigestion()
insilico_digestion_object.setEnzyme(user_enzyme)


def find_all(peptide_seq, protein_seq):
    start = 0
    while True:
        start = protein_seq.find(peptide_seq, start)
        if start == -1: 
            return
        yield start
        start += len(peptide_seq)

def CleaveSiteCheck(peptide_seq, protein_seq, prot_start):
    # Get all occurences of peptide in protein string
    if not prot_start:
        peptide_occurences = list(find_all(peptide_seq, protein_seq))
    else:
        peptide_occurences = [int(prot_start) - 1] #BLAST is 1-based
    ##
    valids = set()
    pep_seq_len = len(peptide_seq)
    ## Check if in proper cleavage site
    for pep_start_index in peptide_occurences:
        if (pep_start_index != -1) and (pep_start_index + pep_seq_len < len(protein_seq)): # Index not before the start of the protein or after start of protein
            ## this function needs as input: the unchanged protein sequence, the position of the peptide overlap, the length of the peptide
            is_valid_cleaved_peptide = insilico_digestion_object.isValidProduct(
                protein_seq,
                pep_start_index,
                pep_seq_len,
                True, # Ignore_missed_cleavages (MC) – Do not compare MC’s of potential peptide to the maximum allowed MC’s
                True, # Allow_nterm_protein_cleavage – Regard peptide as n-terminal of protein if it starts only at pos=1 or 2 (0 or 1 in Python) and protein sequence starts with ‘M’
            )
            valids.add(is_valid_cleaved_peptide)
    ##
    if True in valids:
        is_valid_cleaved_peptide = True
    else:
        is_valid_cleaved_peptide = False
    ##
    return is_valid_cleaved_peptide

protein_seq_match_list = []

known_peptides = set()
counter = 0

outfile_notfound = open(novel_pep_outfile, "w+")
outfile_notfound.write("#Peptide\tNRS\tGENMAP\tGTF\tDarkTrans\tClosestMatchProteinID\tProteinInfo\tMatchInfo\tBestHyperscore\tBestSpectrum\tSecondHyperscore\tSecondSpectrum\n")

outfile_notfound_fasta = open(f"03_EnsembleNotFound_{dataset}.fasta", "w+")

outfile_found = open(f"03_Peptides_Found_Ensembl_{dataset}.txt", "w")
outfile_found.write("#Peptide\tNRS\tGENMAP\tGTF\tDarkTrans\tMatchProteinID\tProteinInfo\tMatchInfo\n")

pep_known_notvalid_cleavage = collections.defaultdict(set)
pep_known_valid_cleavage = collections.defaultdict(set)

# Check for 100% match in database
for pep_seq in pep_set:
    # Regex for Iso/Leucin
    peptide_seq = re.sub("[LI]", "(I|L)", pep_seq)
    # Check the match in database
    for protein_id, protein_seq in ensembl_seqs.items():
        match = re.search(peptide_seq, protein_seq)
        if match:
            known_peptides.add(pep_seq)
            matched_pep_seq = match.group(0)
            # Check cleave site
            valid_cleavage = CleaveSiteCheck(matched_pep_seq, protein_seq, [])
            if valid_cleavage == True:
                pep_known_valid_cleavage[pep_seq].add(protein_id)
            else:
                pep_known_notvalid_cleavage[pep_seq].add(protein_id)

for pep_seq in pep_set:
    MSSamples, HighestHyperscore, SpectrumID, SecondHyperscore, SecondSpectrum, Dark_NRS, AllMatches = PEP_dict[pep_seq][1:]
    if "Thorax" in dataset:
        NRS = "_".join(TRXIDs[Dark_NRS].split("_")[:2])
    else:
        NRS = "_".join(Dark_NRS.split("_")[:2])
    GENMAP = GENMAP_dict[NRS] 
    #
    if pep_seq in known_peptides: #If 100 %
        # And valid cleavage -> found
        if pep_seq in pep_known_valid_cleavage.keys():
            protein_id = ",".join(pep_known_valid_cleavage[pep_seq])
            gene_inf = ",".join([ensembl_ids[prot] for prot in pep_known_valid_cleavage[pep_seq]])
            outfile_found.write("\t".join([pep_seq, NRS, GENMAP, Dark_NRS, protein_id, gene_inf, "Match"]) + "\n")
        else:
            # Not found
            protein_id = ",".join(pep_known_notvalid_cleavage[pep_seq])
            gene_inf = ",".join([ensembl_ids[prot] for prot in pep_known_valid_cleavage[pep_seq]])
            outfile_notfound.write("\t".join([pep_seq, NRS, GENMAP, Dark_NRS, protein_id, gene_inf, "InvalidCleavage", HighestHyperscore, SpectrumID, SecondHyperscore, SecondSpectrum]) + "\n")
    else: # If not 100%, check with blast
        outfile_notfound_fasta.write(f">{pep_seq}\n{pep_seq}\n")

outfile_notfound_fasta.close()

## Blast peptides without 100%
print("Run BLASTp to double check unknown peptides")
#if os.path.isfile(f"03_{dataset}_blastp.m6") == False:
command_line = [f"blast/ncbi-blast-2.14.0+/bin/blastp",
                f"-query 03_EnsembleNotFound_{dataset}.fasta",
                f"-db ensembl_hs_pep/Homo_sapiens.GRCh38.pep.all.fa",
                f"-outfmt '6 qseqid sseqid pident length qstart qend qseq sstart send sseq evalue'",
                f"-out 03_{dataset}_blastp.m6"]
subprocess.call(" ".join(command_line), 
                shell = True)

# Load blast
print("Process BLAST")
blast_matches = collections.defaultdict(list)
low_ims = collections.defaultdict(list)

for line in open(f"03_{dataset}_blastp.m6", "r"):
    qseqid, sseqid, pident, length, qstart, qend, qseq, sstart, send, sseq, evalue = line.split("\t")
    blast_matches[qseqid].append([sseqid, pident, sstart, length, sseq])

for pep, matches in blast_matches.items():
    matches_sorted = sorted(matches, key = lambda x: x[1], reverse = True)
    best_match = matches_sorted[0]
    best_protein_id, highest_pident = best_match[:2]
    ##
    MSSamples, HighestHyperscore, SpectrumID, SecondHyperscore, SecondSpectrum, Dark_NRS, AllMatches = PEP_dict[pep][1:]
    if "Thorax" in dataset:
        NRS = "_".join(TRXIDs[Dark_NRS].split("_")[:2])
    else:
        NRS = "_".join(Dark_NRS.split("_")[:2])
    GENMAP = GENMAP_dict[NRS]
    gene_inf = ensembl_ids[protein_id]
    ##
    if float(highest_pident) > 80:
        # Check if valid cleavage site
        for match in matches_sorted:
            sseqid, pident, sstart, al_len, matched_seq = match
            if float(pident) > 80 and int(al_len)/len(pep) > 0.8: # Just the blast hits >= 80% IM and peptide aligns with more than 80% of it's length
                known_peptides.add(pep)
                protein_sequence = ensembl_seqs[sseqid]
                valid_cleavage = CleaveSiteCheck(pep, protein_sequence, sstart)
                if valid_cleavage == True:
                    if valid_cleavage == True:
                        pep_known_valid_cleavage[pep_seq].add(sseqid)
                    else:
                        pep_known_notvalid_cleavage[pep_seq].add(sseqid)
    ##
    if pep in known_peptides:
        if pep in pep_known_valid_cleavage:
            protein_id = ",".join(pep_known_valid_cleavage[pep_seq])
            gene_inf = ",".join([ensembl_ids[prot] for prot in pep_known_valid_cleavage[pep_seq]])
            outfile_found.write("\t".join([pep, NRS, GENMAP, Dark_NRS, protein_id, gene_inf, f"BLAST_Match"]) + "\n")
        elif pep in pep_known_notvalid_cleavage:
            protein_id = ",".join(pep_known_notvalid_cleavage[pep])
            gene_inf = ",".join([ensembl_ids[prot] for prot in pep_known_valid_cleavage[pep]])
            outfile_notfound.write("\t".join([pep, NRS, GENMAP, Dark_NRS, protein_id, gene_inf, f"InvalidCleavage,IM:{highest_pident}", HighestHyperscore, SpectrumID, SecondHyperscore, SecondSpectrum]) + "\n")
    else:
        gene_inf = ensembl_ids[best_protein_id]
        outfile_notfound.write("\t".join([pep, NRS, GENMAP, Dark_NRS, best_protein_id, gene_inf, f"IM:{highest_pident}", HighestHyperscore, SpectrumID, SecondHyperscore, SecondSpectrum]) + "\n")
