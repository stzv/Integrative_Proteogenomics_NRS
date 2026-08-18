============================================================
README — NRS RNA-SEQ ANALYSIS PIPELINE
============================================================

SHORT DESCRIPTION
------------------------------------------------------------
This pipeline evaluates transcriptional support for
non-reference sequences (NRS) using publicly available
RNA-seq datasets.

The workflow processes RNA-seq alignments against an
extended reference genome containing both GRCh38 and NRS,
identifies NRS supported by RNA-seq reads, computes
population-level expression frequencies, estimates
expression levels using counts per million (CPM), and
integrates transcriptional support with genomic placement
information.

The resulting datasets provide:

- RNA-seq support for NRS
- Per-sample NRS coverage information
- Population frequency of NRS expression
- Expression levels (CPM)
- Tissue-specific expression patterns
- Genomic annotation of expressed NRS
- Integration of expression and genomic placement data

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

Python packages

  - pandas
  - numpy
  - matplotlib
  - seaborn
  - pysam
  - more_itertools
  - BioPython
  - collections
  - regex

External tools

  - STAR
  - Samtools

Reference files

  - Extended reference genome (GRCh38 + NRS)
  - NRS FASTA assembly
  - RNA-seq datasets
  - SDRF metadata files
  - NRS genomic placement information

------------------------------------------------------------
INPUT DATA
------------------------------------------------------------

This repository contains analysis scripts only.

The original input files are not included because they are
stored in restricted institutional storage.

Required inputs include:

  - Public RNA-seq datasets
  - STAR alignment results
  - NRS FASTA assembly
  - NRS coverage tables
  - Sample metadata (SDRF files)
  - NRS genomic location information
  - Gene annotations

File names used in the scripts should be treated as
examples and may need adjustment for local installations.

------------------------------------------------------------
FOLDER LAYOUT
------------------------------------------------------------

NRS_RNASeq/

  00_extract_rna_domain.py
  01_merge_LOG_results.py
  02_depth_analysis.py
  03_saturation_curve_randomized.py
  04_process_coverage.py
  05_merge_CJ.py
  06_CJ_annotation.py
  07_NRS_RNA_table_v2.py
  07a_tissue_expression_picture.py
  07b_BoxPlot.py
  08_Check_Reads_per_domain_v2.py
  09_RNASeq_GenMapping.py

------------------------------------------------------------
PIPELINE OVERVIEW
------------------------------------------------------------

(1) Extract RNA-specific NRS domains

    Script:
      00_extract_rna_domain.py

    Purpose:
      Extract NRS sequence domains relevant for RNA-seq
      expression analyses.

    Output:
      RNA-domain sequence files.


(2) Collect RNA-seq mapping statistics

    Script:
      01_merge_LOG_results.py

    Purpose:
      Merge STAR alignment statistics across RNA-seq
      samples.

      Collect numbers of uniquely mapped reads used for
      downstream CPM normalization.

    Output:
      RNA-seq sample statistics table.


(3) Determine RNA-seq support for NRS

    Script:
      02_depth_analysis.py

    Purpose:
      Calculate read depth across NRS for each RNA-seq
      sample.

      Generate per-sample NRS support information based on
      continuous read coverage.

    Output:
      02_SABE_1172_UNHESMSV_NRS_RNA_reads_count.txt


(4) Generate saturation curves

    Script:
      03_saturation_curve_randomized.py

    Purpose:
      Estimate cumulative discovery of expressed NRS with
      increasing numbers of RNA-seq samples.

    Output:
      Saturation curve tables and figures.


(5) Calculate expression frequencies

    Script:
      04_process_coverage.py

    Purpose:
      Convert per-sample support into NRS expression
      frequencies across the complete RNA-seq dataset.

    Output:
      RNA-seq frequency tables.


(6) Integrate genomic mappings

    Script:
      05_merge_CJ.py

    Purpose:
      Merge RNA-supported NRS with genomic placement
      information.

    Output:
      Combined location-expression tables.


(7) Annotate expressed NRS

    Script:
      06_CJ_annotation.py

    Purpose:
      Annotate expressed NRS relative to genes and genomic
      features.

    Output:
      Annotated expression tables.


(8) Calculate CPM values and tissue-specific expression

    Script:
      07_NRS_RNA_table_v2.py

    Purpose:
      Compute counts per million (CPM) for all expressed
      NRS and integrate RNA sample metadata.

      Generate tissue-level and sample-level expression
      summaries.

    Output:

      07_SABE_1172_UNHESMSV_NRS_RNA_tissue_table.txt

      07_SABE_1172_UNHESMSV_NRS_RNA_tissue_table_cpm>10.txt

      07_SABE_1172_UNHESMSV_PublicRNA_selection.txt


(9) Summarize tissue-level expression

    Script:
      07a_tissue_expression_picture.py

    Purpose:
      Compute average tissue expression and generate
      tissue-specific expression summaries.

    Output:
      Tissue expression plots.


(10) Compare expression between tissues

    Script:
      07b_BoxPlot.py

    Purpose:
      Generate tissue-level expression statistics and
      visualization tables.

    Output:
      Tissue expression summary tables.


(11) Evaluate RNA-seq read support

    Script:
      08_Check_Reads_per_domain_v2.py

    Purpose:
      Assess RNA-seq read support across RNA-supported
      domains and evaluate expression robustness.

    Output:
      Read-support summaries.


(12) Integrate expression and GENMAP locations

    Script:
      09_RNASeq_GenMapping.py

    Purpose:
      Combine RNA-seq expression support with final NRS
      genomic location assignments.

    Output:
      Final expression-location integration tables.

------------------------------------------------------------
DECISION RULES
------------------------------------------------------------

RNA-seq coverage support

  Mapping quality:

      MAPQ >= 20

  Read support:

      >= 10 reads

  Continuity requirement:

      >= 10 consecutive bp

  Reads containing >=10 bp soft clipping at either end
  are removed before coverage analysis.

------------------------------------------------------------

Expression frequency

  An NRS is considered expressed in a sample when it
  satisfies the RNA-seq coverage criteria above.

  Expression frequency is calculated as the proportion of
  RNA-seq samples supporting the NRS.

------------------------------------------------------------

CPM calculation

  CPM values are calculated using the number of uniquely
  mapped reads reported by STAR.

  For each NRS, the maximum CPM observed across all
  samples is retained for downstream summaries.

------------------------------------------------------------

Tissue expression

  Tissue identity is extracted from SDRF metadata.

  Expression summaries are generated separately for each
  tissue and tissue category.

------------------------------------------------------------
OUTPUT FILES
------------------------------------------------------------

Main outputs:

  02_SABE_1172_UNHESMSV_NRS_RNA_reads_count.txt

  RNA-seq saturation curve files

  RNA-seq expression frequency tables

  Annotated RNA expression tables

  07_SABE_1172_UNHESMSV_NRS_RNA_tissue_table.txt

  07_SABE_1172_UNHESMSV_NRS_RNA_tissue_table_cpm>10.txt

  Tissue-expression figures

  Expression-GENMAP integration tables

------------------------------------------------------------
RUN ORDER
------------------------------------------------------------

  00_extract_rna_domain.py

  01_merge_LOG_results.py

  02_depth_analysis.py

  03_saturation_curve_randomized.py

  04_process_coverage.py

  05_merge_CJ.py

  06_CJ_annotation.py

  07_NRS_RNA_table_v2.py

  07a_tissue_expression_picture.py

  07b_BoxPlot.py

  08_Check_Reads_per_domain_v2.py

  09_RNASeq_GenMapping.py

------------------------------------------------------------
NOTES
------------------------------------------------------------

- This pipeline uses publicly available RNA-seq datasets
  to evaluate transcriptional support for NRS.

- RNA-seq alignments are performed against an extended
  reference genome containing both GRCh38 and NRS.

- Input data are not included in this repository.

- Several files originate from restricted institutional
  storage and cannot be redistributed publicly.

- Absolute file paths present in the original scripts
  should be replaced by local paths before execution.

============================================================
END OF README
============================================================
