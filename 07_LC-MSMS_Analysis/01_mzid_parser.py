from glob import glob
import collections
import pastaq
from glob import glob
import os

import sys

mzid_files = glob("SwisProt/SwissProt_ORF_k81/*_SwissProt_DPD_normalSearchMSFragger_DP*_report.mzid")
all_spectrum_files = glob(f"SwisProt/SwissProt_ORF_k81/*_SwissProt_DPD_normalSearchMSFragger_DP*_psm.tsv")

header = list()
accession_dict = collections.defaultdict(dict)
spectrum_dict = collections.defaultdict(dict)
peptides = collections.defaultdict(dict)

for fil in mzid_files:
  header = list()
  accession_dict = collections.defaultdict(dict)
  # Get sample ID
  sampleid = fil.split("/")[-1].split("_")[0]
  print("Processing", sampleid)
  # Parse mzid file
  mzid_object = pastaq.read_mzidentml(str(fil)) # 'db_sequences', 'peptide_evidence', 'peptides', 'spectrum_matches']
  # Fetch sequences information
  db_sequence_info = mzid_object.db_sequences
  # Keep trinity sequences
  for entry in db_sequence_info:
    if str(entry.accession).startswith("TRINITY"):
      accession_dict[entry.accession] = set()
  # Fetch peptide evidence
  bd_evidence = mzid_object.peptide_evidence
  for pep in bd_evidence:
    # Filter out decoy and non-trinity sequences
    if pep.decoy == 0 and pep.peptide_id and pep.db_sequence_id.startswith("TRINITY"):
      accession_dict[pep.db_sequence_id].add(pep.peptide_id)
  # Fetch spectrum matches
  spec_fil = f"SwisProt/SwissProt_ORF/{sampleid.lower()}_SwissProt_DPD_normalSearchMSFragger_DP{sampleid.lower()}_psm.tsv"
  for fil in all_spectrum_files: # To do case insensitive search for the file => not uniform capitals in naming of files
    if spec_fil.lower() == fil.lower():
      spectrum_file = fil
  for line in open(spectrum_file, "r"):
    # Skip header
    if line.startswith("Spectrum"):
      continue
    # 150210_2_01_c.04116.04116.4	Control2_SwissProt_DPD_normalSearchMSFragger_DPcontrol2.pep.xml	HHEEEIVHHKK
    #
    spectrum_id, fil_loc, peptide = line.split()[:3]
    score = line.split()[17]
    if peptide not in spectrum_dict:
      spectrum_dict[peptide] = set()
    spectrum_dict[peptide].add(f"{spectrum_id}:{score}:{sampleid}")
  # Count how many trinity seqs have peptide evidence
  supported_count = [evid for evid in accession_dict.values() if evid]
  header.append(f"{sampleid}:{len(supported_count)}")
  #
  outfile = open(f"01_supported_ORF_{sampleid}.txt", "w+")
  for seq, evid in accession_dict.items():
    if evid:
      outfile.write(f"{seq}\t{','.join(evid)}\n")
      for p in evid:
        if not "Samples" in peptides[p].keys():
          peptides[p]["Samples"] = set()
        if not "Seqs" in peptides[p].keys():
          peptides[p]["Seqs"] = set()
        if not "MSIDs" in peptides[p].keys():
          peptides[p]["MSIDs"] = ""
        ##
        peptides[p]["Samples"].add(sampleid)
        peptides[p]["Seqs"].add(seq)
        peptides[p]["MSIDs"] = spectrum_dict[p]


pepout = open("01_peptides_list.txt", "w+")
pepout.write("#Peptide\tSamplesCount\tSeqCount\tSamples\tSeqs\tMSSpectraCount\tMSSpectraIDs\n")
for p in peptides.keys():    
  pepout.write(f"{p}\t{len(peptides[p]['Samples'])}\t{len(peptides[p]['Seqs'])}\t{','.join(peptides[p]['Samples'])}\t{','.join(peptides[p]['Seqs'])}\t{len(peptides[p]['MSIDs'])}\t{','.join(peptides[p]['MSIDs'])}\n")
