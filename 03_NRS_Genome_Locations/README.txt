============================================================
README — GENOMIC PLACEMENT OF NON-REFERENCE SEQUENCES
(GENMAP PIPELINE)
============================================================

SHORT DESCRIPTION
------------------------------------------------------------
This pipeline assigns genomic locations to non-reference
sequences (NRS) assembled from discordant sequencing reads.

The workflow integrates four independent mapping strategies:

- Partial Mapping Coordinates (PMC) derived from Minimap2
  alignments
- Anchoring using paired-end read links between NRS and
  GRCh38
- 10X Linked Reads using barcode-coherent long-range
  information
- Chimeric RNA-Seq alignments linking NRS to reference
  genomic locations

These mapping signals are integrated using a rule-based
decision framework to generate a single consensus genomic
location per NRS, referred to as the GENMAP coordinate.

The resulting datasets provide:

- Genomic placement evidence from multiple independent
  mapping approaches
- Consensus genomic locations for NRS
- Genomic interval assignments standardized across
  mapping methods
- Functional annotation of mapped NRS relative to genes,
  exons, introns, and intergenic regions

Content of README:
- Project information
- Author
- Dependencies
- Input data information
- Folder layout
- Pipeline overview
- Mapping methods
- GENMAP decision framework
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
  - pandas >= 1.3
  - numpy >= 1.21
  - matplotlib >= 3.4
  - matplotlib-venn >= 0.11
  - upsetplot >= 0.6

External tools
  - Minimap2
  - Bowtie2
  - STAR
  - Samtools

Reference files
  - GRCh38 reference genome
  - Ensembl / GENCODE GTF annotation
  - Final NRS FASTA assembly

------------------------------------------------------------
INPUT DATA
------------------------------------------------------------

This repository contains analysis scripts only.

The original input files are not included because they are
stored in restricted institutional storage.

Required inputs include:

  - Final NRS FASTA assembly
  - Minimap2 alignments against GRCh38
  - 10X Linked Read mapping results
  - Anchor read mappings
  - Chimeric RNA-Seq mappings
  - GRCh38 reference genome
  - Ensembl / GENCODE gene annotation

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

02_Consensus_GENMAP_incl_ChimericReads/

  01_merge_locations_into_table.py
  01b_compare_locsv4.py
  02_compare_CJ_mappingv4.py
  03_GENMAP_vs_GTF.py

------------------------------------------------------------
PIPELINE OVERVIEW
------------------------------------------------------------

PART I — GENERATION OF LOCATION EVIDENCE

The first part of the workflow generates candidate genomic
locations using three independent mapping approaches.

(1) 10X Linked Read processing

    Scripts:
      01_merging_10x_locs.py
      02_count_10x_reads.py

    Purpose:
      - Merge genomic locations detected across all 10X
        samples
      - Cluster nearby locations
      - Count supporting reads, barcodes, and samples
      - Select the strongest-supported 10X location


(2) Anchoring analysis

    Scripts:
      03a_get_anchors.pl
      03b_count_anchor_reads.py
      07_matched_best_anchor.py

    Purpose:
      - Detect paired-end links between NRS and GRCh38
      - Group nearby anchor positions
      - Calculate read support
      - Determine insertion orientation
      - Select the most robust anchor location


(3) Partial Mapping Coordinates (PMC)

    Script:
      06_PMC_supporting_10X_or_anchor.py

    Purpose:
      - Incorporate evidence from partially mapped NRS
      - Compare PMC locations to 10X and Anchoring
        evidence
      - Identify concordant placements


(4) Comparison of mapping methods

    Scripts:
      05_10x_vs_anchors.py
      05b_matched_loc_lengths.py
      08_Robust_10x_Anchor_finalloc.py

    Purpose:
      - Compare locations inferred by different methods
      - Detect concordant placements
      - Generate a robust placement table for downstream
        integration

PART II — CONSENSUS GENMAP ASSIGNMENT

The second part of the workflow combines all available
mapping evidence, including Chimeric RNA-Seq mappings.

(5) Merge all location evidence

    Script:
      01_merge_locations_into_table.py

    Purpose:
      Create a unified table containing:

      - 10X placements
      - Anchoring placements
      - PMC placements
      - Chimeric RNA-Seq placements


(6) Compare locations between methods

    Script:
      01b_compare_locsv4.py

    Purpose:
      Detect concordant location pairs between mapping
      methods using predefined distance thresholds.


(7) Generate consensus GENMAP assignment

    Script:
      02_compare_CJ_mappingv4.py

    Purpose:
      Apply GENMAP decision rules and select a single
      genomic position per NRS.


(8) Annotate genomic locations

    Script:
      03_GENMAP_vs_GTF.py

    Purpose:
      Compare final GENMAP locations to gene models and
      determine overlap with:

      - Genes
      - Exons
      - Introns
      - Intergenic regions

------------------------------------------------------------
MAPPING METHODS
------------------------------------------------------------

1. Partial Mapping Coordinates (PMC)

PMC provides high-resolution evidence but is sensitive to
repetitive sequence.

Logic:

  - NRS with >=95% identity and >=200 bp inserted or
    soft-clipped segments are classified as partially
    mapped.

Stored information:

  - Chromosome
  - Coordinate
  - MAPQ
  - CIGAR string
  - Alignment metadata

------------------------------------------------------------

2. Anchoring

Anchoring provides high positional precision.

Logic:

  - Anchor positions are grouped into clusters within
    <=500 bp.
  - Clusters supported by >=5 read pairs are retained.
  - Directional consistency is used to select the best
    anchor.

Final precision:

  Approximately +/-500 bp.

------------------------------------------------------------

3. 10X Linked Reads

10X placements provide lower positional precision but
strong genome-wide support.

Logic:

  - Nearby placements are clustered.
  - Read, barcode, and sample support are summarized.
  - The strongest-supported cluster is selected.

Final precision:

  Approximately +/-50 kb.

------------------------------------------------------------

4. Chimeric RNA-Seq

Available only for expressed NRS.

Logic:

  - Chimeric NRS↔GRCh38 links are retained unless they
    originate from the HLA region.
  - Breakpoint midpoint is used as the coordinate.
  - RNA placements are compared to other methods using
    a 50 kb concordance threshold.

------------------------------------------------------------
GENMAP DECISION FRAMEWORK
------------------------------------------------------------

STEP 1 — Detect matching location pairs

A matching pair is defined as two independent mapping
methods placing the same NRS in approximately the same
genomic region.

Distance thresholds:

  10X ↔ Anchoring       <= 50 kb
  10X ↔ PMC             <= 50 kb
  10X ↔ RNA-Seq         <= 50 kb
  Anchoring ↔ PMC       <= 1 kb
  Anchoring ↔ RNA-Seq   <= 50 kb
  PMC ↔ RNA-Seq         <= 50 kb

------------------------------------------------------------

STEP 2 — If exactly one matching pair exists

  10X + Anchoring
      → Anchoring

  10X + RNA-Seq
      → RNA-Seq

  Anchoring + PMC
      → Anchoring

  PMC + RNA-Seq
      → RNA-Seq

------------------------------------------------------------

STEP 3 — If zero or multiple matching pairs exist

Fallback hierarchy:

  1. Use 10X location if available
  2. Else use a unique Anchoring location
  3. Else mark NRS as unmapped

------------------------------------------------------------

STEP 4 — Convert placement into GENMAP interval

Anchoring

  - +/-500 bp around anchor position
  - Format:
      chr:start:end:strand

Chimeric RNA-Seq

  - Breakpoint-derived coordinates
  - Format:
      chr:start:end:*

10X Linked Reads

  - +/-50 kb around linked-read position
  - Format:
      chr:start:end:*

Mitochondrial, decoy, and viral sequences

  - Retained unchanged
  - Excluded from downstream analyses

------------------------------------------------------------
OUTPUT FILES
------------------------------------------------------------

Main outputs:

  Robust 10X / Anchoring placement table

  Integrated location evidence table

  Matched location pairs

  Consensus GENMAP assignment table

  GENMAP gene annotation table

  Location summary statistics

------------------------------------------------------------
RUN ORDER
------------------------------------------------------------

PART I — Location Evidence Generation

  01_merging_10x_locs.py

  02_count_10x_reads.py

  03a_get_anchors.pl

  03b_count_anchor_reads.py

  05_10x_vs_anchors.py

  05b_matched_loc_lengths.py

  06_PMC_supporting_10X_or_anchor.py

  07_matched_best_anchor.py

  08_Robust_10x_Anchor_finalloc.py

PART II — Consensus Assignment

  01_merge_locations_into_table.py

  01b_compare_locsv4.py

  02_compare_CJ_mappingv4.py

  03_GENMAP_vs_GTF.py

------------------------------------------------------------
NOTES
------------------------------------------------------------

- Individual NRS can receive support from multiple
  mapping methods.

- GENMAP produces a single final genomic placement per
  NRS.

- Chimeric RNA-Seq mappings are only available for
  transcriptionally active NRS.

- Input data are not included in this repository.

- Several files originate from restricted institutional
  storage and cannot be redistributed publicly.

- Absolute file paths present in the original scripts
  should be replaced by local paths before execution.

============================================================
END OF README
============================================================
