============================================================
README — NRS ASSEMBLY AND PARTIAL MAPPING ANALYSIS PIPELINE
============================================================

SHORT DESCRIPTION
------------------------------------------------------------
This pipeline was developed to characterize non-reference
sequences (NRS) assembled from discordant sequencing reads
and evaluate their relationship to the GRCh38 reference
genome.

The workflow identifies fully unmapped and partially mapped
NRS, generates a high-confidence final NRS ("dark freeze")
set, summarizes population-level coverage across samples,
and detects candidate chimeric NRS through analysis of
primary and supplementary alignments.

The resulting datasets provide:
- A final collection of non-reference sequences
- Genomic placement evidence for partially mapped NRS
- Per-sample NRS coverage information
- NRS population frequencies
- Candidate chimeric contigs requiring additional
  interpretation

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
Integrative proteogenomic analysis of non-reference sequences

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

Languages
  - Python 3.8+
  - Perl

Python packages
  - Biopython
  - numpy
  - regex
  - collections

External tools
  - Minimap2
  - Samtools
  - GNU core utilities
    (grep, awk, zcat)

------------------------------------------------------------
INPUT DATA
------------------------------------------------------------

This repository contains analysis scripts only.

The original input files are not included because they are
stored in restricted institutional storage.

Required inputs include:
  - NRS FASTA assembly
  - GRCh38 reference genome
  - Minimap2 alignments
  - Contig length tables
  - Coverage and frequency tables generated during earlier
    processing steps

File names used in the scripts should be treated as examples
and may need adjustment for local installations.

------------------------------------------------------------
FOLDER LAYOUT
------------------------------------------------------------

01_NRS_Assembly/

  03_minimap2_dark_vs_GRCh38.pl
  04_split_alignment.py
  05_make_partial_map.py
  06_extract_NRS_freeze.py
  07_extract_NRS_freeze_coverage.py
  08_listNRS_coverage.py
  09_NRS_chimeric_connections.py

------------------------------------------------------------
PIPELINE OVERVIEW
------------------------------------------------------------

(1) Align NRS contigs against GRCh38
    Script:
      03_minimap2_dark_vs_GRCh38.pl
    Purpose:
      Align assembled NRS contigs to the GRCh38 reference
      genome using Minimap2.
    Input:
      - NRS FASTA
      - GRCh38 reference genome
    Output:
      SABE_1172_UNHESMSV_vs_GRCh38.sam.gz

(2) Classify contigs by alignment status
    Script:
      04_split_alignment.py
    Purpose:
      Separate contigs into:
        - mapped
        - partially mapped
        - unmapped
      based on alignment identity and size 
      of the non-reference segment.

(3) Generate partial mapping coordinates
    Script:
      05_make_partial_map.py
    Purpose:
      Extract genomic placement evidence from partially
      mapped contigs.
    Output:
      Partial Mapping Coordinates (PMC) for downstream
      location analyses.

(4) Generate final NRS ("dark freeze") set
    Script:
      06_extract_NRS_freeze.py
    Purpose:
      Create the final high-confidence NRS collection for
      downstream analyses.
    Output:
      SABE_1172_UNHESMSV_NRS_dark_freeze_final.fa

(5) Generate NRS coverage matrix
    Script:
      07_extract_NRS_freeze_coverage.py
    Purpose:
      Combine sample-level coverage information for all
      NRS in the final freeze set.
    Output:
      SABE1172_UNHESMSV_NRS_dark_freeze_final_coverage.txt

(6) Calculate NRS population frequencies
    Script:
      08_listNRS_coverage.py
    Purpose:
      Determine the fraction of samples supporting each
      NRS.
    Output:
      SABE_1172_UNHESMSV_NRSList_IndividList.txt

(7) Resolve candidate chimeric NRS
    Script:
      09_NRS_chimeric_connections.py
    Purpose:
      - Identify partially mapped NRS containing
      supplementary alignments and determine whether the
      supplementary alignments explain previously
      unmapped sequence.
      - Contigs whose supplementary alignments collectively
      resolve the non-reference region are flagged
      separately from contigs that retain substantial
      unmapped sequence.
    Output:
      09_NRS_SecMapping_resolved.txt

------------------------------------------------------------
DECISION RULES
------------------------------------------------------------

Partial mapping classification

  A contig is considered partially mapped when:
    - sequence divergence is less than 5%
    - and at least 200 bp remains as an inserted or
      unmatched segment

Unmapped classification

  A contig is considered unmapped when:
    - sequence divergence is 5% or greater

Chimeric alignment evaluation

  For partially mapped contigs:

    - primary and supplementary alignments are combined
    - aligned regions are reconstructed using CIGAR
      operations
    - contigs are retained as unresolved if >= 200 bp
      remains unmapped after considering supplementary
      alignments

------------------------------------------------------------
OUTPUT FILES
------------------------------------------------------------

Main outputs:

  SABE_1172_UNHESMSV_NRS_dark_freeze_final.fa
  SABE1172_UNHESMSV_NRS_dark_freeze_final_coverage.txt
  SABE_1172_UNHESMSV_NRSList_IndividList.txt
  09_NRS_SecMapping_resolved.txt

------------------------------------------------------------
RUN ORDER
------------------------------------------------------------

  03_minimap2_dark_vs_GRCh38.pl
  04_split_alignment.py
  05_make_partial_map.py
  06_extract_NRS_freeze.py
  07_extract_NRS_freeze_coverage.py
  08_listNRS_coverage.py
  09_NRS_chimeric_connections.py

------------------------------------------------------------
NOTES
------------------------------------------------------------
- Input data are not included in this repository.
- Several files originate from restricted institutional
  storage and cannot be redistributed publicly.
- Absolute file paths present in the original scripts
  should be replaced by local paths before execution.

============================================================
END OF README
============================================================
