#!/usr/bin/perl -w
use strict;

system join( ' ', 'minimap2', '-ax', 'sr',
                              '-t', 24, 
                              'Brazil/ref_1172/blast_GRCh38/GRCh38.p13.alt.hs38d1.fa', 
                              'megahit_UNHESMSV_merged/SABE1172_unmapped.fa', 
                  '|', 'gzip -c >SABE_1172_UNHESMSV_vs_GRCh38.sam.gz',
           );
