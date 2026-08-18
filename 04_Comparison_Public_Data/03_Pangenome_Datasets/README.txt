============================================================
README — PANGENOME DATASET COMPARISON PIPELINE
============================================================

SHORT DESCRIPTION
------------------------------------------------------------
This pipeline compares SABE non-reference sequences (NRS)
against publicly available human pangenome datasets.

The workflow aligns NRS contigs to external assemblies,
classifies the resulting alignments, and extracts matching
contigs and genomic positions for downstream comparison.

The resulting datasets provide:

- Identification of NRS already represented in external
  pangenome resources
- Genomic placements of matching NRS
- Classification of complete and partial matches
- Comparison between SABE NRS and public pangenome
  assemblies

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

Languages
  - Python 3.8+
  - Perl

Python packages
  - gzip
  - glob

External tools
  - Minimap2
  - Samtools

Reference files
  - Public pangenome assemblies
  - Final NRS FASTA assembly

------------------------------------------------------------
INPUT DATA
------------------------------------------------------------

This repository contains analysis scripts only.

The original input files are not included because they are
stored in restricted institutional storage.

Required inputs include:

  - Final NRS FASTA assembly
  - Public pangenome reference assemblies
  - Minimap2 alignment files

File names used in the scripts should be treated as
examples and may need adjustment for local installations.

------------------------------------------------------------
FOLDER LAYOUT
------------------------------------------------------------

03_Pangenome_Datasets/

  03_minimap2_dark_vs_publicdata.pl
  04_split_sr_w_overlaps.pl
  09_extract_seqnames.py
  10_extract_seqnames_blast.py

------------------------------------------------------------
PIPELINE OVERVIEW
------------------------------------------------------------

(1) Align NRS against public pangenome datasets

    Script:
      03_minimap2_dark_vs_publicdata.pl

    Purpose:
      Align NRS contigs to external pangenome assemblies
      using Minimap2.

    Output:
      SAM alignment files.


(2) Classify alignment status

    Script:
      04_split_sr_w_overlaps.pl

    Purpose:
      Categorize alignments into mapped and partially
      mapped groups according to alignment overlap and
      completeness.

    Output:
      Filtered alignment files.


(3) Extract aligned NRS and coordinates

    Script:
      09_extract_seqnames.py

    Purpose:
      Extract NRS identifiers, lengths, chromosomes,
      and genomic coordinates from mapped and partially
      mapped alignment files.

    Output:
      Lists of aligned NRS with genomic coordinates.

(4) Extract matching contig identifiers

    Script:
      10_extract_seqnames_blast.py

    Purpose:
      Summarize sequence matches identified in comparisons
      against public pangenome assemblies.

    Output:
      Lists of matching NRS and reference contigs.

------------------------------------------------------------
DECISION RULES
------------------------------------------------------------

Alignment classification

  - NRS are classified according to mapping status
    after alignment to public pangenome datasets.

  - Complete and partial mappings are retained
    separately.

Uniqueness filtering

  - Duplicate NRS identifiers are removed from the
    final summary tables.

------------------------------------------------------------
OUTPUT FILES
------------------------------------------------------------

Main outputs:

  Mapped NRS tables

  Partially mapped NRS tables

  NRS coordinate lists

  Matched public pangenome contig lists

------------------------------------------------------------
RUN ORDER
------------------------------------------------------------

  03_minimap2_dark_vs_publicdata.pl

  04_split_sr_w_overlaps.pl

  09_extract_seqnames.py

  10_extract_seqnames_blast.py

------------------------------------------------------------
NOTES
------------------------------------------------------------

- This pipeline evaluates the overlap between SABE NRS
  and publicly available pangenome resources.

- Input data are not included in this repository.

- Several files originate from restricted institutional
  storage and cannot be redistributed publicly.

- Absolute file paths present in the original scripts
  should be replaced by local paths before execution.

============================================================
END OF README
============================================================
