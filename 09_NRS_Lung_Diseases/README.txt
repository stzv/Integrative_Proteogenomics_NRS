============================================================
README — NRS LUNG DISEASE ASSOCIATION PIPELINE
============================================================

SHORT DESCRIPTION
------------------------------------------------------------
This pipeline investigates the relationship between
non-reference sequence (NRS) expression and chronic lung
disease across multiple respiratory cohorts.

The workflow evaluates associations between NRS expression
and disease status using case-control analyses, computes
odds ratios for NRS occurrence in disease versus control
samples, and performs differential expression analysis
using edgeR. Significant NRS and annotated genes are then
visualized using volcano plots, heatmaps, and clustering
analyses.

The resulting datasets provide:

- Cohort-specific NRS disease associations
- Odds ratio estimates for NRS expression
- Differential expression results for genes and NRS
- Disease-associated NRS candidates
- Heatmap-based expression clusters
- Cross-cohort comparison of disease expression trends

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
  - R

Python packages

  - pandas
  - numpy
  - matplotlib
  - seaborn
  - upsetplot

R packages

  - edgeR
  - pheatmap
  - dplyr
  - ggplot2
  - ggrepel
  - biomaRt
  - RColorBrewer

------------------------------------------------------------
INPUT DATA
------------------------------------------------------------

This repository contains analysis scripts only.

The original input files are not included because they are
stored in restricted institutional storage.

Required inputs include:

  - NRS expression frequency tables
  - RNA-seq count matrices
  - Patient phenotype information
  - Cohort metadata
  - Ensembl gene annotations

File names used in the scripts should be treated as
examples and may need adjustment for local installations.

------------------------------------------------------------
FOLDER LAYOUT
------------------------------------------------------------

09_NRS_Lung_Diseases/

  01_NRS_vs_Individuals_v3.py

  02a_Differential_Expression1.Rmd

  02b_Differential_Expression_graphs.py

------------------------------------------------------------
PIPELINE OVERVIEW
------------------------------------------------------------

(1) Identify disease-associated NRS

    Script:
      01_NRS_vs_Individuals_v3.py

    Purpose:

      Integrate patient metadata and NRS expression
      information across cohorts.

      For each NRS:

      - Determine presence across patients
      - Calculate population frequency
      - Calculate disease/control counts
      - Calculate odds ratios
      - Compare disease and control enrichment

    Output:

      Patient information tables

      NRS disease association tables

      Odds-ratio distributions

      Cohort comparison plots


(2) Differential expression analysis

    Script:
      02a_Differential_Expression1.Rmd

    Purpose:

      Perform differential expression analysis using
      edgeR.

      Disease and control groups are compared while
      correcting for:

      - Sex
      - Age

      Differential expression is computed for both:

      - Annotated genes
      - Non-reference sequences (NRS)

    Output:

      Differential expression tables

      Significant DE genes

      Significant DE NRS


(3) Differential expression visualization

    Script:
      02b_Differential_Expression_graphs.py

    Purpose:

      Generate summary plots and visualizations from
      differential expression analyses.

    Output:

      Publication-ready figures and expression summaries

------------------------------------------------------------
DECISION RULES
------------------------------------------------------------

Disease association analysis

  An NRS must occur in at least 10% of patients within
  a cohort to be included in odds-ratio analyses.

  Odds ratios are calculated using disease and control
  patient counts for each NRS.

  Continuity correction (+0.5) is applied to avoid
  division by zero.

------------------------------------------------------------

Differential expression filtering

  A feature is retained if:

    CPM >= 1

    in at least half of one of the compared groups

  OR

    NRS total read count >= 100

------------------------------------------------------------

Differential expression model

  edgeR quasi-likelihood framework

  Model:

    Expression ~ Disease Status + Sex + Age

------------------------------------------------------------

Significance criteria

  Differential expression significance:

    FDR <= 0.05

------------------------------------------------------------

Gene annotation

  Significant genes are annotated using Ensembl /
  BioMart and linked to:

    - HGNC symbols
    - Gene descriptions

------------------------------------------------------------
OUTPUT FILES
------------------------------------------------------------

Main outputs:

  Patient_Information_All_Cohorts.txt

  NRS disease association tables

  Odds-ratio summary figures

  Differential expression tables

  Significant DE genes

  Significant DE NRS

  Volcano plots

  Heatmaps

  Expression cluster tables

  Clustered NRS summaries

------------------------------------------------------------
RUN ORDER
------------------------------------------------------------

  01_NRS_vs_Individuals_v3.py

  02a_Differential_Expression1.Rmd

  02b_Differential_Expression_graphs.py

------------------------------------------------------------
NOTES
------------------------------------------------------------

- This pipeline evaluates associations between NRS
  expression and lung disease phenotypes.

- Disease enrichment analyses are performed separately
  for each cohort.

- Differential expression analyses simultaneously
  evaluate annotated genes and NRS.

- Input data are not included in this repository.

- Several files originate from restricted institutional
  storage and cannot be redistributed publicly.

- Absolute file paths present in the original scripts
  should be replaced by local paths before execution.

============================================================
END OF README
============================================================
