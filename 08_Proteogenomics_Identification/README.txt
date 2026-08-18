============================================================
README — PUBLIC COHORT PROTEOGENOMIC SCREEN (ARMS, PRESTO, THORAX)
============================================================

This folder contains scripts used to identify predicted ORF-derived
peptides in public LC–MS/MS datasets from ARMS, Presto, and Thorax
cohorts. The workflow includes MSFragger searches, peptide evidence
aggregation, novelty assessment against Ensembl, peptide-to-ORF mapping,
NRS positional analysis, cohort merging, and visualization.


This folder contains scripts used to identify predicted NRS‑derived
ORF peptides in large publicly available proteomic datasets.
The workflow builds directly on the ORF‑prediction and validation
pipeline. The same predicted ORFs and TransDecoder‑derived peptide sequences 
generated from the RNA-Seq asthma and COPD cohorts (referred to as ARMS, Presto and Thorax) 
were used here to screen two large ProteomeXchange datasets for evidence of
translation of non‑reference genomic sequence.

Content of README:
- Project information
- Author
- Dependencies
- Input data information
- Folder layout
- Pipeline overview
- Decision rules
- Output files
- Run order

------------------------------------------------------------
PROJECT INFORMATION
------------------------------------------------------------

Research project:

  Integrative proteogenomic analysis of non-reference
  sequences

Related publication:

  DOI: TO BE ADDED AFTER PUBLICATION

------------------------------------------------------------
AUTHOR / CONTACT
------------------------------------------------------------

Stepanka Zverinova
University Medical Center Groningen (UMCG)
ORCID: 0000-0002-3370-9484

Peter Horvatovich
Groningen University
0000-0003-2218-1140

------------------------------------------------------------
DEPENDENCIES
------------------------------------------------------------

Python >= 3.8
  glob
  collections
  BioPython (SeqIO)
  gzip
  re
  subprocess
  pandas

R
  ggplot2
  MSnbase
  mzR
  protViz
  stringr

External Tools
  MSFragger
  Philosopher
  BLASTp
  SwissProt protein database
  Ensembl GRCh38 peptide database

------------------------------------------------------------
PIPELINE OVERVIEW
------------------------------------------------------------

(1) Download public MS datasets  
    Script: 00_download_data.py
    - Downloads raw MS files from PXD020192 and PXD026370.

(2) Generate input and sample metadata  
    Script: 00a_Get_Sample_Names.py
    - Extracts mzML filenames and sample identifiers for MSFragger batch runs.

(3) Run MSFragger + Philosopher  
    Script: 01_MSFraggerPhilosopherAnalysis_Thorax.py
    - Searches mzML files against SwissProt + ORF-derived peptides.
    - Produces mzIdentML and PSM tables.
    - Pipeline reused for ARMS and Presto cohorts.

(4) Extract ORF-supported (“dark”) peptides  
    Script: 02_Get_Dark_peptides_v2_<cohort>.py
    - Loads PSM tables.
    - Retains peptides mapping to NRS or TRX ORFs.
    - Keeps highest-scoring PSM per peptide.
    - Output: 02_peptides_list_<cohort>.txt

(5) Novelty assessment against Ensembl (exact match + BLASTp)
    Script: 03_Check_Novel_Peptides_v4.py
    - Exact match search in Ensembl + SwissProt (I/L ambiguity allowed).
    - BLASTp for remaining peptides.
    - Novel peptide criteria:
         * <80% identity, AND
         * <80% coverage
    - Outputs:
         03_Peptides_Found_Ensembl_<cohort>.txt
         03_Peptides_NotFound_Ensembl_<cohort>.txt
         03_<cohort>_blastp.m6

(6) Peptide → ORF → NRS positional mapping  
    Scripts: 03a_Peptide_Position_<cohort>.py and 06_Pep_Location.py
    - Maps peptide positions onto ORFs.
    - Computes NRS-derived regions from ORF→NRS alignments.
    - Identifies peptides located within NRS-derived ORF sequence.

(7) Merge ARMS, Presto, and Thorax results  
    Script: 04_Merge_Cohorts.py
    - Combines novel peptides across cohorts.
    - Merges NRS IDs, GENMAP hits, sample counts, PSM scores.
    - Outputs:
         04_NovelPeptides_Merged.tsv
         04_Novel_peptides_list_publicdata.txt
         04_NovelPeptides.bed

(8) Peptide visualization  
    Script: 05_Run_Visualization.py
    - Calls associated R scripts to plot annotated MS/MS spectra for
      selected peptides.

------------------------------------------------------------
DECISION RULES
------------------------------------------------------------

- Keep only PSMs mapping to NRS ORF sequences.
- Retain highest Hyperscore per peptide per cohort.
- Peptide is KNOWN if:
      * exact match found in Ensembl/SwissProt, OR
      * BLASTp match ≥80% identity AND ≥80% coverage.
- Peptide is NOVEL if neither condition above is met.
- Peptide considered NRS-derived if peptide coordinates lie within the
  NRS-overlapping region of the ORF.

------------------------------------------------------------
OUTPUT FILES
------------------------------------------------------------

02_peptides_list_<cohort>.txt  
03_Peptides_NotFound_Ensembl_<cohort>.txt  
03_Peptides_Found_Ensembl_<cohort>.txt  
03_<cohort>_blastp.m6  
03a_Peptide_Position_on_ORF.tsv  
04_NovelPeptides_Merged.tsv  
04_NovelPeptides.bed  

------------------------------------------------------------
RUN ORDER
------------------------------------------------------------

00_download_data.py  
00a_Get_Sample_Names.py  
01_MSFraggerPhilosopherAnalysis_Thorax.py  
02_Get_Dark_peptides_v2_<cohort>.py  
03_Check_Novel_Peptides_v4.py  
03a_Peptide_Position_<cohort>.py  
04_Merge_Cohorts.py  
05_Run_Visualization.py  
06_Pep_Location.py

============================================================
END OF README
============================================================
