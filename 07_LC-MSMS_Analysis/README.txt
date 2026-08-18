============================================================
README — LC‑MS/MS ANALYSIS OF NOVEL ORF‑DERIVED PEPTIDES

This folder documents the workflow used to identify and validate
proteomic evidence for predicted ORFs containing non‑reference segments.
The pipeline includes MS database searching, PSM extraction, 
aggregation of supporting evidence per ORF, filtering against the 
Ensembl human proteome, positional mapping of peptides onto ORFs and
their NRS-derived regions, and visualization of endogenous and synthetic
MS/MS spectra.

CONTENT
  DEPENDENCIES
  FOLDER LAYOUT
  PIPELINE OVERVIEW
  DECISION RULES (SUMMARY)
  OUTPUT FILES
  RUN ORDER
============================================================

------------------------------------------------------------
DEPENDENCIES
------------------------------------------------------------

Python 3.8 or higher
  glob
  gzip
  collections
  BioPython (SeqIO)
  re
  numpy
  pandas
  matplotlib

R
  mzR
  MSnbase
  protViz
  ggplot2
  stringr

External Tools
  MSFragger (search engine)
  mzIdentML (.mzid) output
  mzML raw files
  SwissProt protein database
  Custom ORF peptide FASTA (from TransDecoder)

------------------------------------------------------------
FOLDER LAYOUT
------------------------------------------------------------

Main analysis scripts:

  01_mzid_parser.py
  02_merge_evidence.py
  03_search_peptide_ensembl.py
  05_Find_Peptide_Position.py
  06_MSMS_Visualization_SZ.R
  07_plot_annotated_spectrum.R
  07a_Spectra_functions.R

Intermediate files include peptide lists, ORF support tables,
Ensembl novelty results, NRS positional mapping tables, and annotated
spectral plots.

------------------------------------------------------------
PIPELINE OVERVIEW
------------------------------------------------------------

(1) Parse MSFragger mzIdentML results and collect peptide evidence
    Script: 01_mzid_parser.py
    - Parses .mzid files for peptide-spectrum matches (PSMs).
    - Retains only non‑decoy PSMs that map to TRINITY-derived ORFs.
    - Extracts peptide IDs, sample IDs, sequence accessions, and
      corresponding spectrum identifiers and scores.
    - Outputs:
         01_supported_ORF_<sample>.txt
         01_peptides_list.txt

(2) Aggregate ORF-level evidence across all samples
    Script: 02_merge_evidence.py
    - Loads grouped ORFs (TRX IDs).
    - For each ORF group:
         * number of samples with peptide evidence
         * number of unique peptides
         * total supporting spectra
         * MS-support frequency (samples / 19)
    - Outputs:
         02_Merged_Peptide_Evidence.txt
         02_Predicted_NotSupported.txt
         02_List_All_Supported_ORFs.fa
    - Produces summary plots of population vs proteomic support.

(3) Compare observed peptides against Ensembl human peptides
    Script: 03_search_peptide_ensembl.py
    - Loads Ensembl GRCh38 peptide database.
    - Searches for exact matches (I/L interchangeable).
    - Peptides found in Ensembl are classified as known.
    - Others are retained as candidate novel peptides.
    - Outputs:
         03_Peptides_Found_Ensembl.txt
         03_Peptides_NotFound_Ensembl.txt

(4) Map peptide locations onto ORFs and NRS-derived regions
    Script: 05_Find_Peptide_Position.py
    - Loads TRX ORF sequences and ORF-to-NRS alignment metadata.
    - Reconstructs NRS-overlapping segments using CIGAR operations.
    - Determines peptide positions on protein coordinates.
    - Evaluates whether each peptide overlaps the NRS-derived
      (non-reference) region in the ORF.
    - Outputs:
         05_UnknownPeptides_Position_ORF.txt

(5) Visualize MS/MS spectra for endogenous and synthetic peptides
    Scripts: 06_MSMS_Visualization_SZ.R,
             07_plot_annotated_spectrum.R,
             07a_Spectra_functions.R
    - Loads mzML files for endogenous and synthetic peptides.
    - Extracts MS/MS spectra for best-scoring PSM per peptide.
    - Predicts fragment ion series (b- and y-ions).
    - Generates annotated mirror plots to compare endogenous and
      synthetic fragmentation patterns.
    - Used for visual validation of peptide identity.

------------------------------------------------------------
DECISION RULES (SUMMARY)
------------------------------------------------------------

(1) PSM filtering
      - Only non‑decoy PSMs are accepted.
      - Only PSMs assigned to TRINITY-derived ORFs are retained.
      - For peptides with multiple PSMs in a sample, the highest scoring
        match is used.

(2) ORF support criteria
      - ORF considered proteomically supported if ≥1 peptide is observed.
      - Support frequency is calculated from 19 peptide-MS samples.

(3) Peptide novelty criteria
      - Peptide must not appear in the Ensembl GRCh38 peptide database.
      - Peptide location must overlap the NRS-derived region of its ORF
        (protein coordinate mapping).

(4) ORF novelty criteria
      - ORF passes earlier filtering (mapping quality, divergence, etc.).
      - ORF does not align to GRCh38 extended reference with high identity
        (tested before LC-MS/MS analysis).

(5) Spectral confirmation
      - Endogenous MS/MS spectra must closely match synthetic peptide
        spectra in mirror plots.

------------------------------------------------------------
OUTPUT FILES
------------------------------------------------------------

  01_supported_ORF_<sample>.txt
  01_peptides_list.txt
  02_Merged_Peptide_Evidence.txt
  02_Predicted_NotSupported.txt
  02_List_All_Supported_ORFs.fa
  03_Peptides_Found_Ensembl.txt
  03_Peptides_NotFound_Ensembl.txt
  05_UnknownPeptides_Position_ORF.txt
  Mirror-plot spectra files (PDF/JPEG)

------------------------------------------------------------
RUN ORDER
------------------------------------------------------------

  01_mzid_parser.py
  02_merge_evidence.py
  03_search_peptide_ensembl.py
  05_Find_Peptide_Position.py
  06_MSMS_Visualization_SZ.R
  07_plot_annotated_spectrum.R

============================================================
END OF README
============================================================
