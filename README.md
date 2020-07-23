# Non-Reference Sequences (NRS)
## Introduction


***************************************************
## Stage I - Assembly

## Stage II - NRS selection
System requirements: **Include links to the packages**
- Perl
- Bio::SeqIO
- Bowtie2
- Minimap2

How it works:
*01_rename_merge_index.pl*
Contigs resulting from the Stage I Assembly are tagged as NRS and merged with the reference genome.
Bowtie2 creates indexed files.

*02_bowtie2_vs_extended.pl*
FastQ files containing the paired-end reads of targed genomes are aligned to the extended reference (reference genome + NRS).
The amount of reads mapping to NRS with specified mapping quality (default = 20) are saved per sample in *bowtie_vs_extended/{sample}_dark_cnt.txt.gz*.
Anchoring information is saved per sample in *{sample}_anchors.txt.gz*.

*03_minimap2_dark_vs_GRCh38.pl*
Tagged NRS is aligned to reference genome.

*04_split_sr_w_overlaps.pl* **Need script for creating contig_len.txt.gz file**
The NRS aligned to reference genome are split into i) unmapped NRS, ii) mapped NRS, and iii) partially mapped RNS. 
- Unmapped NRS: Identity match < 95%
- Mapped NRS: Identity match >= 95%
- Partially mapped NRS: Identity match >= 95% & length of unmapped segment >= 200bp

Identity match and unmapped segment length variables can be adjusted.
The file 04_map.txt.gz contains the binary map of the matches and insertions between NRS and reference, based on the CIGAR string. **Please check**

*05_make_partial_map.pl*
Check for presence/absence of unmapped and partially mapped NRS in the samples. Keeps NRS with coverage >= 7.5X & =< 100X, calculated as ***reads count * reads length / NRS length***. 
Frequency of NRS in samples is saved in *05_freq_results.txt*

*06_extract_dark_freeze.pl*
Based on binary map, mapped segments of partially mapped NRS are cut out. Where applicable, the mapping information of the mapped segments is stored.
The unmapped NRS and unmapped segments of partially mapped NRS are saved in *NRS_dark_freeze_coverage.txt*

## Stage III - Genome mapping
