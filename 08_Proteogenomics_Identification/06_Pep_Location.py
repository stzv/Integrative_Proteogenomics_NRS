import re
import collections
import math

import sys

## Load in ORF sequences
TRXids = collections.defaultdict(dict)

for line in open("00_Predicted_ORFs_grouped_IDs.txt", "r"):
 if line.startswith("#"):
  continue
 #
 TRX_ID,	PopFreq,	All_IDs,	SEQ, PepSeq = line.rstrip().split("\t")
 ids = [id.split(":")[1] for id in All_IDs.split(",")]
 TRXids[TRX_ID] = PepSeq

## Load spectrums per peptide
spectrums = collections.defaultdict(dict)
for line in open("01_peptides_list.txt", "r"):
   pep = line.rstrip().split("\t")[0]
   ids  = line.rstrip().split("\t")[-1]
   spectrums[pep] = ids

## Load in ORF alignment
NRS = dict()

for line in open("00_ORF_Alignment.txt", "r"):
 nrs, nrs_len, trxid, alignment = line.rstrip().split("\t")
 sample, orfid, tl, fl, pos, cg, div, al_seq, og_seq = alignment.split("|")
 # Get NRS loc coord on peptide seq
 # Add values
 value = {"nrs_id": nrs, "trinity_seq_len": tl.replace("tl:", ""), "fl": fl.replace("fl:", ""), "alignment_start": int(pos.replace("pos:", "")) - 1, "cigar": cg.replace("cg:", ""), "pred_pep_seq": og_seq, "peptide_seq": TRXids[trxid]}
 NRS.setdefault(trxid, list())
 if value not in NRS[trxid]:
  NRS[trxid].append(value)

## Load in unknown peptides
NotFoundPeptides = collections.defaultdict(dict)

for line in open("03_Peptides_NotFound_Ensembl.txt", "r"):
 # Skip header
 if line.startswith("#"):
  continue
 ##
 Peptide, PredictedSeqFreq, EvidenceFrequence, TRXIds, Orfs = line.rstrip().split("\t")
 # Find the positions
 for trxid in TRXIds.split(","):
  orfseq = TRXids[trxid]
  pep_loc = f"{orfseq.find(Peptide)}:{int(orfseq.find(Peptide)) + len(Peptide) - 1}"
  addition = {"Peptide": Peptide, "TRXId": trxid, "PepLoc": pep_loc, "NRSLoc": "", "NRSLoc_Pep": ""}
  #
  NotFoundPeptides.setdefault(trxid, list())
  if addition not in NotFoundPeptides[trxid]:
    NotFoundPeptides[trxid].append(addition)

#######
outfile = open("05_UnknownPeptides_Position_ORF.txt", "w+")
outfile.write("#Peptide\tPepLen\tOrfId\tOrfLen\tPepWithinNRS\tPepOnOrfPos\tNrsOnOrfPosNt\tNrsOnOrfPosPep\tORFSeq\tMSSpectrums\n")

for k, v in NotFoundPeptides.items():
  for entry in v:
    for e in NRS[k]:
      cg_split = re.findall("\d+\w", e["cigar"])
      if int(e["fl"]) == 16:
        cg_split = cg_split[::-1]
      nrs_overlap = ""
      for cg in cg_split:
        if cg[-1] == "S":
          nrs_overlap = nrs_overlap + "0" * int(cg[:-1])
        elif cg[-1] in ["M", "I"]:
          nrs_overlap = nrs_overlap + "1" * int(cg[:-1])
        elif cg[-1] == "D": # Deletion -> NRS extra information that is missing in reference => Skip
          continue
      nrs_start = nrs_overlap.find("1")
      nrs_end = nrs_overlap.rfind("1")
      nrs_loc = f"{nrs_start}:{nrs_end}"
      nrs_start_pep = int(nrs_start/3)
      nrs_end_pep = int(nrs_end/3)
      ## Cannot get completely precise nt to pep position due to possible framshift => int rounding allows for 1 aa error
      nrs_loc_pep = f"{nrs_start_pep}:{nrs_end_pep}"
      ## Find if NRS overlaps in it's aa location the position of the evidence peptide
      pep_location = entry["PepLoc"].split(":")
      mx = max(nrs_start_pep, int(pep_location[0]))
      mn = min(nrs_end_pep, int(pep_location[1]))
      overlap = range(mx, mn)
      if overlap:
        pep_within_nrs = "Yes"
      else:
        pep_within_nrs = "No"
      ## Note NRS part of ORF in small letters
      orf = TRXids[k]
      orf = orf.replace(orf[nrs_start_pep:nrs_end_pep], orf[nrs_start_pep:nrs_end_pep].lower())
      #
      if pep_within_nrs == "Yes":
        outfile.write(f'{entry["Peptide"]}\t{len(entry["Peptide"])}\t{k}\t{len(TRXids[k])}\t{pep_within_nrs}\t{entry["PepLoc"]}\t{nrs_loc}\t{nrs_loc_pep}\t{e["fl"]}\t{orf}\t{spectrums[entry["Peptide"]]}\n')
      #
