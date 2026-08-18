============================================================
README — NRS COMPARISON AGAINST PUBLIC DATASETS AND CHM13
============================================================

SHORT DESCRIPTION
------------------------------------------------------------
This pipeline compares non-reference sequences (NRS) from
the SABE cohort against publicly available genomic
resources, including human pangenome datasets, RefSeq
non-redundant sequences, and the CHM13v2+Y reference
assembly.

The workflow identifies NRS already represented in public
resources, evaluates overlap between datasets, and
determines which NRS can only be resolved using CHM13.

The resulting datasets provide:

- Identification of NRS already represented in public
  pangenome resources
- Identification of NRS matching RefSeq non-reference
  sequences
- Identification of NRS uniquely resolved by CHM13
- Comparative analysis of sequence sharing across
  population datasets
- Population frequency characterization of shared and
  SABE-specific NRS

Content of README:
- Project information
- Author
- Pipeline overview
- Decision-making summary
- Input files
- Comparison workflows
- Folder layout
- Dependencies
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
PIPELINE OVERVIEW
------------------------------------------------------------

This workflow compares SABE NRS against several external
resources to determine whether the detected sequences have
been previously observed in other studies or reference
assemblies.

Datasets evaluated:

  - HUPAN
  - PDBU01
  - GoNL
  - RefSeq Non-Redundant (NR)
  - CHM13v2+Y

Pipeline summary:

  1. Align SABE NRS against each public dataset using
     Minimap2.

  2. Classify alignments into:
       - Mapped
       - Partially mapped
       - Unmapped

  3. Extract NRS detected in each dataset.

  4. Compare NRS against RefSeq NR using BLAST.

  5. Integrate evidence across datasets.

  6. Identify:
       - Shared NRS
       - Dataset-specific NRS
       - SABE-unique NRS
       - CHM13-only placements

  7. Visualize overlap between datasets using UpSet
     analyses and summary statistics.

------------------------------------------------------------
DECISION-MAKING SUMMARY
------------------------------------------------------------

Alignment classification

  Unmapped:
      divergence (de:f tag) > 0.05

  Partially mapped:
      CIGAR contains a segment >= 200 bp composed of
      insertion, deletion, or soft-clipped sequence

  Mapped:
      all remaining primary alignments

Presence in a public dataset

  An NRS is considered present if it is:
      - mapped
      - partially mapped

  Unmapped contigs are considered absent.

RefSeq NR criteria

  percent identity >= 99%

  alignment coverage >= 95% of NRS length

CHM13-only NRS

  NRS appears in the CHM13 mapped list

  AND

  NRS does NOT appear in any GENMAP category:
      - PMC
      - Anchoring
      - 10X
      - RNA-Seq

  AND

  NRS is not classified as completely unlocalized.

  Therefore, genomic placement support originates
  exclusively from CHM13.

No additional filtering

  - No coverage filters
  - No realignment
  - No polishing
  - No identity threshold beyond alignment
    classification criteria

------------------------------------------------------------
INPUT FILES
------------------------------------------------------------

Main inputs

  SABE NRS FASTA

      SABE1172_UNHESMSV_NRS_dark_freeze_final.fa

  Frequency table

      SABE_1172_UNHESMSV_frequency_nrs.txt

Public datasets

  HUPAN
      DOI: 10.1186/s13059-019-1751-y

  PDBU01
      DOI: 10.1038/s41588-018-0273-y

  GoNL
      DOI: 10.1038/ejhg.2013.118

RefSeq NR

  Input:
      SABE_1172_UNHESMSV_NCBINR.m6

CHM13

  Input:
      03_Chm13v2_SABE_UNHESMSV_mapped_list.txt

------------------------------------------------------------
COMPARISON WORKFLOWS
------------------------------------------------------------

01_Chm13

  Purpose:

    Determine whether SABE NRS absent from GRCh38 can be
    localized using the complete CHM13v2+Y assembly.

  Main question:

    Are some SABE NRS derived from genomic regions that
    are missing or unresolved in GRCh38 but present in
    CHM13?

------------------------------------------------------------

02_NCBI_nonredundant

  Purpose:

    Compare SABE NRS against the NCBI non-redundant
    nucleotide database.

  Main question:

    Have these sequences already been reported in public
    sequence repositories?

------------------------------------------------------------

03_Pangenome_Datasets

  Purpose:

    Compare SABE NRS against publicly available human
    pangenome resources.

  Main question:

    Are these sequences shared across human populations
    and represented in existing pangenome assemblies?

------------------------------------------------------------

Together, these comparisons allow distinction between:

  - NRS already represented in public databases
  - NRS present in other population datasets
  - NRS resolved only by CHM13
  - Potentially novel SABE-specific sequences

Further details for each comparison workflow are provided
in the README files within the corresponding pipeline
folders.

------------------------------------------------------------
FOLDER LAYOUT
------------------------------------------------------------

01_Chm13/

    00_minimap2.py
    01_prep_files.py
    02_split_samfile.pl
    03_Extract_Seqnames.py
    04_compare_CJ_Chm13Y_mappingv3.py
    NRS_only_Chm13.Rmd

------------------------------------------------------------

02_NCBI_nonredundant/

    blast_vs_nr.py

------------------------------------------------------------

03_Pangenome_Datasets/

    03_minimap2_dark_vs_publicdata.pl
    04_split_sr_w_overlaps.pl
    09_extract_seqnames.py
    10_extract_seqnames_blast.py
    D_upsetplot_v2.py

------------------------------------------------------------
DEPENDENCIES
------------------------------------------------------------

Python (3.8+)

    pandas
    numpy
    matplotlib
    upsetplot
    glob
    gzip

Perl

    No non-standard Perl modules required.

R

    karyoploteR
    GenomicRanges
    dplyr
    tidyr
    stringr
    scales
    grid
    RColorBrewer
    corrplot
    circlize

External tools

    Minimap2
        Used with -ax sr preset for alignment of SABE
        NRS against public datasets.

    BLAST
        Used for comparison against RefSeq NR.

    gzip
        Used for compressed SAM processing.

------------------------------------------------------------
NOTES
------------------------------------------------------------

- This repository contains analysis scripts only.

- Public datasets referenced in the workflow are
  distributed by their original authors and are not
  included in this repository.

- Several intermediate files originate from restricted
  institutional storage and cannot be redistributed.

- Absolute file paths present in the original scripts
  should be replaced by local paths before execution.

- The NCBI non-redundant comparison workflow currently
  contains the BLAST search step only and may be expanded
  when additional processing scripts become available.

============================================================
END OF README
============================================================
