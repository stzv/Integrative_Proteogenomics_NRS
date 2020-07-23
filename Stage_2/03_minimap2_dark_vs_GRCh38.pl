#!/usr/bin/perl -w
use strict;

## ENVIRONMENT VARIABLES ##
my $hg_reference = 'GRCh38_full_analysis_set_plus_decoy_hla.fa'
my $threads = 24


# Alignment NRS to Reference
system join( ' ', './minimap2', '-ax', 'sr',
                              '-t', $threads, 
                              $hg_reference, 
                              'NRS.fa', 
                  '|', 'gzip -c >NRS_vs_GRCh38.sam.gz',
           );
