============================================================
README — ORF Prediction and Proteomic Confirmation
============================================================

This folder documents the workflow used to identify novel open reading
frames (ORFs) within NRS-containing transcripts and to confirm their
translation using proteomic evidence. The pipeline consists of transcript
assembly, ORF prediction, NRS‑alignment filtering, ORF grouping,
reference‑genome comparison, and preparation of a peptide search library
for LC‑MS/MS matching.

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
- Notes

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

------------------------------------------------------------
DEPENDENCIES
------------------------------------------------------------

Python >= 3.8
  glob
  subprocess
  BioPython (SeqIO)
  gzip
  collections
  numpy
  matplotlib

External tools:
  Trinity (v2.13.2 or similar)
  TransDecoder (LongOrfs, Predict)
  minimap2 (alignments of ORFs to NRS)
  BLAST+ (megablast)
  Tar/Gzip (packaging peptide FASTAs)

------------------------------------------------------------
FOLDER ORGANIZATION
------------------------------------------------------------

This folder contains all scripts used for transcript assembly,
ORF prediction, ORF→NRS alignment filtering, ORF grouping,
reference-genome comparison, and ID renaming. The main scripts are:

  01_transcript_assembly.py
  02_predict_ORF.py
  03_transcript_alignment.py
  04_transcript_alignment_process.py
  05_separate_fasta.py
  07_Group_same_ORFs.py
  08_ORF_alignment_rename_TRXID.py
  09_BLAST_ORF_GRCh38.py

Intermediate output files include transcript assemblies, predicted ORFs,
alignment files, grouped ORFs, and final TRX‑annotated ORF tables.

------------------------------------------------------------
PIPELINE OVERVIEW
------------------------------------------------------------

(1) De novo transcript assembly
    Script: 01_transcript_assembly.py  
    - Runs Trinity separately for each RNA‑Seq sample.
    - Outputs: *_trinity.Trinity.fasta

(2) ORF prediction from transcript assemblies
    Script: 02_predict_ORF.py
    - Runs TransDecoder.LongOrfs and TransDecoder.Predict.
    - Generates predicted ORF nucleotide sequences (*.cds)
      and translated peptide sequences (*.pep).

(3) Alignment of ORFs to the NRS catalog
    Script: 03_transcript_alignment.py
    - Aligns predicted ORF CDS sequences to the SABE NRS reference
      using minimap2 (-ax sr).
    - Outputs compressed SAM files for each sample.

(4) Filtering ORFs based on alignment quality
    Script: 04_transcript_alignment_process.py 
    Applied criteria:
        - Mapping quality <20 → discard
        - Supplementary alignments → discard
        - Sequence divergence >2% → discard
    Retained ORFs are saved with alignment details, including:
        transcript ID, CIGAR, divergence, aligned sequence, and
        original ORF sequence.

(5) Extracting ORF peptides for each individual
    Script: 05_separate_fasta.py
    - Looks up predicted peptide sequences for ORFs passing all filters.
    - Writes one FASTA per individual (Predicted_ORF_inclNRS/<ID>_ORF.fasta).
    - Produces a tar.gz archive of all ORF peptide FASTAs.

(6) Grouping ORFs across individuals
    Script: 07_Group_same_ORFs.py
    - ORFs are grouped by identical CDS sequence.
    - Computes population frequency of each ORF group.
    - Creates:
        07_PredictedORFs_grouped.fa           (peptides)
        07_PredictedORFs_grouped_cds.fa       (CDS)
        07_PredictedORFs_grouped_IDs.txt      (ID list)
    - Generates population‑frequency histograms.

(7) Assigning TRX IDs to ORF alignments
    Script: 08_ORF_alignment_rename_TRXID.py
    - Maps each ORF alignment entry to its corresponding TRX group.
    - Produces: 08_THORAX_SABE1172UNHESMSV_ORF_alignment_TRXIDs.txt

(8) Assessing ORF novelty relative to GRCh38 extended reference
    Script: 09_BLAST_ORF_GRCh38.py
    - Runs BLASTN (megablast) of grouped ORF CDS sequences against
      the extended GRCh38 reference.
    - ORFs with high‑identity BLAST matches (≥95%) are flagged as
      likely reference‑derived and can be removed from downstream analyses.

------------------------------------------------------------
DECISION RULES (SUMMARY)
------------------------------------------------------------

(1) ORF alignment acceptance criteria:
      - Mapping quality ≥20
      - No supplementary alignments
      - Divergence ≤2% (identity ≥98%)
      All others are discarded.

(2) ORF grouping:
      - ORFs sharing identical CDS sequences are grouped together.
      - Population frequency = (# individuals containing ORF) / total

(3) Novelty filtering:
      - BLASTN identity <95% to GRCh38 extended reference
        → retain as potentially non-reference ORF. 

(4) Peptide library generation:
      - Peptides from retained ORFs (grouped) form the basis of the
        custom proteogenomic library for LC‑MS/MS searches.

------------------------------------------------------------
OUTPUTS
------------------------------------------------------------

Primary outputs include:

  *_trinity.Trinity.fasta
  *_trinity.Trinity.fasta.transdecoder.cds
  *_trinity.Trinity.fasta.transdecoder.pep
  *_ORF_SABE1172UNHESMSV_alignment.sam.gz
  04_THORAX_SABE1172UNHESMSV_ORF_alignment.txt
  Predicted_ORF_inclNRS/<ID>_ORF.fasta
  Predicted_ORF_inclNRS.tar.gz
  07_PredictedORFs_grouped.fa
  07_PredictedORFs_grouped_cds.fa
  07_PredictedORFs_grouped_IDs.txt
  08_THORAX_SABE1172UNHESMSV_ORF_alignment_TRXIDs.txt
  09_Thorax_ORF_GRh38_extended_IM95.txt

These files serve as the transcript‑level and peptide‑level inputs for
the proteogenomic search library used in LC‑MS/MS analysis.

============================================================
END OF README
============================================================
