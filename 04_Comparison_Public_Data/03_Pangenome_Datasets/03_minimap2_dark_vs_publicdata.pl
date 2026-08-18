#!/usr/bin/perl -w
use strict;

#my $setname = $ARGV[0]; #'PDBU01'; #'GoNL'; #'HAN';

#die if $#ARGV == -1;

print "Running Public database overlap";

my $fastafile = "SABE_1172_UNHESMSV_genotyping/SABE1172_UNHESMSV_NRS_dark_freeze_final.fa";

print " HAN";
my $setname = "HAN";
system join( ' ', 'minimap2', '-ax', 'sr', '-t', 24, 
                              'SABE_1172_UN_complete/ref_1172/public_dark/'.$setname.'.fa', 
                              $fastafile, 
                  '|', 'gzip -c > 03_SABE1172_UNHESMSV_vs_'.$setname.'.sam.gz');

print " PDBU01";
my $setname = "PDBU01";

system join( ' ', 'minimap2', '-ax', 'sr', '-t', 24, 
		'SABE_1172_UN_complete/ref_1172/public_dark/'.$setname.'.fa',
		$fastafile,
		'|', 'gzip -c > 03_SABE1172_UNHESMSV_vs_'.$setname.'.sam.gz');

print " GoNL";
my $setname = "GoNL";

system join( ' ', 'minimap2', '-ax', 'sr', '-t', 24,
		'SABE_1172_UN_complete/ref_1172/public_dark/'.$setname.'.fa',
		$fastafile,
		'|', 'gzip -c > 03_SABE1172_UNHESMSV_vs_'.$setname.'.sam.gz');
