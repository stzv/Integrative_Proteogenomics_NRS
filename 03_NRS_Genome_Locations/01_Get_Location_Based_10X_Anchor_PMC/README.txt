============================================================
README — GENMAP LOCATION EVIDENCE GENERATION PIPELINE
============================================================

SHORT DESCRIPTION
------------------------------------------------------------
This pipeline generates genomic placement evidence for
non-reference sequences (NRS) using three independent
mapping strategies:

- 10X Linked Reads
- Anchoring (paired-end read linking)
- Partial Mapping Coordinates (PMC)

The workflow integrates genomic signals from linked-read
barcodes, paired-end reads, and partial alignments to
identify candidate genomic locations for NRS. Each method
is optimized according to its strengths and precision,
resulting in a robust set of candidate placements for
downstream integration in the GENMAP consensus framework.

The resulting datasets provide:
- Candidate genomic locations based on 10X Linked Reads
- Candidate genomic locations based on Anchoring
- Supporting evidence from partial mappings (PMC)
- Concordance between mapping strategies
- Final high-confidence genomic placements supported by
  10X, Anchoring, and PMC evidence

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
  - pandas
  - numpy
  - matplotlib
  - collections

External tools
  - Bowtie2
  - Minimap2
  - Samtools
  - gzip

------------------------------------------------------------
INPUT DATA
------------------------------------------------------------

This repository contains analysis scripts only.

The original input files are not included because they are
stored in restricted institutional storage.

Required inputs include:

  - Final NRS FASTA assembly
  - 10X Linked Read mapping results
  - Bowtie2 anchor mappings
  - Partial Mapping Coordinate (PMC) tables
  - Contig length information
  - GRCh38 genomic coordinates

File names used in the scripts should be treated as
examples and may need adjustment for local installations.

------------------------------------------------------------
FOLDER LAYOUT
------------------------------------------------------------

01_GENMAP_Based_10X_Anchor_PMC/

  01_merging_10x_locs.py
  02_count_10x_reads.py
  03a_get_anchors.pl
  03b_count_anchor_reads.py
  05_10x_vs_anchors.py
  05b_matched_loc_lengths.py
  06_PMC_supporting_10X_or_anchor.py
  07_matched_best_anchor.py
  08_Robust_10x_Anchor_finalloc.py

------------------------------------------------------------
PIPELINE OVERVIEW
------------------------------------------------------------

(1) Merge 10X Linked Read locations

    Script:
      01_merging_10x_locs.py

    Purpose:
      Merge NRS genomic locations detected across all
      10X Linked Read samples into a unified table.

    Output:
      Combined 10X location evidence for each NRS.


(2) Summarize 10X support

    Script:
      02_count_10x_reads.py

    Purpose:
      Cluster nearby 10X locations and calculate:

      - Supporting reads
      - Supporting barcodes
      - Number of supporting samples

      The highest-supported genomic location is retained
      as the primary 10X placement.

    Output:
      Ranked 10X candidate locations.


(3) Generate Anchoring evidence

    Script:
      03a_get_anchors.pl

    Purpose:
      Extract paired-end anchor signals connecting NRS
      and GRCh38 locations.

      Anchor direction and insertion orientation are
      determined from paired-end mapping information.

    Output:
      Raw anchor evidence table.


(4) Summarize Anchoring support

    Script:
      03b_count_anchor_reads.py

    Purpose:
      Group anchor positions into location clusters and
      count supporting read pairs.

    Output:
      Candidate anchor locations with support counts.


(5) Compare 10X and Anchoring placements

    Script:
      05_10x_vs_anchors.py

    Purpose:
      Identify NRS for which 10X and Anchoring predict
      concordant genomic locations.

      Matching placements are identified using predefined
      genomic distance thresholds.

    Output:
      Matched and unmatched location pairs.


(6) Evaluate location support lengths

    Script:
      05b_matched_loc_lengths.py

    Purpose:
      Summarize contig characteristics and evaluate the
      size distribution of matched mapping events.

    Output:
      Supporting descriptive statistics.


(7) Add Partial Mapping Coordinate evidence

    Script:
      06_PMC_supporting_10X_or_anchor.py

    Purpose:
      Compare PMC locations with 10X and Anchoring
      evidence and identify concordant placements.

    Output:
      PMC-supported mapping evidence.


(8) Select best Anchoring placement

    Script:
      07_matched_best_anchor.py

    Purpose:
      Resolve multiple anchor candidates and select the
      most robust anchor location based on strand and
      directional consistency.

    Output:
      Best-anchor assignments.


(9) Generate final robust placements

    Script:
      08_Robust_10x_Anchor_finalloc.py

    Purpose:
      Integrate 10X, Anchoring, and PMC evidence into a
      final set of high-confidence genomic locations.

    Output:
      Robust genomic placement table.

------------------------------------------------------------
DECISION RULES
------------------------------------------------------------

10X Linked Reads

  - Nearby locations are clustered.
  - Supporting reads, barcodes, and sample counts are
    summarized.
  - The strongest-supported cluster is retained as the
    primary 10X location.

Anchoring

  - Anchor positions are clustered within 500 bp.
  - Clusters supported by at least five read pairs are
    retained.
  - Directional consistency is used to resolve multiple
    anchor candidates.

10X ↔ Anchoring Concordance

  - Locations are considered concordant when they occur
    within 50 kb of each other.

Partial Mapping Coordinates (PMC)

  - PMC evidence is used as supporting evidence for
    confirming candidate locations.
  - NRS must satisfy predefined partial-mapping
    criteria established during NRS assembly.

Robust Placements

  - Concordant locations supported by multiple mapping
    strategies are prioritized over single-method
    placements.

------------------------------------------------------------
OUTPUT FILES
------------------------------------------------------------

Main outputs:

  Merged 10X location table

  10X support summary table

  Anchor support summary table

  Matched 10X-anchor locations

  PMC-supported location table

  Best-anchor assignments

  Final robust genomic location table

------------------------------------------------------------
RUN ORDER
------------------------------------------------------------

  01_merging_10x_locs.py

  02_count_10x_reads.py

  03a_get_anchors.pl

  03b_count_anchor_reads.py

  05_10x_vs_anchors.py

  05b_matched_loc_lengths.py

  06_PMC_supporting_10X_or_anchor.py

  07_matched_best_anchor.py

  08_Robust_10x_Anchor_finalloc.py

------------------------------------------------------------
NOTES
------------------------------------------------------------

- This pipeline generates placement evidence only.

- Final location assignment and integration of Chimeric
  RNA-Seq mappings are performed in the downstream
  GENMAP consensus pipeline.

- Input data are not included in this repository.

- Several files originate from restricted institutional
  storage and cannot be redistributed publicly.

- Absolute file paths present in the original scripts
  should be replaced by local paths before execution.

============================================================
END OF README
============================================================
