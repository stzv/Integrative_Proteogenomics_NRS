#!/usr/bin/perl -w
use strict;
use lib '/usr/local/ensembl72/bioperl-live';
use Bio::SeqIO;

### ENVIRONMENT VARIABLES - change as necessary
## Tagging dark matter DNA
my $contigs_in = 'contigs.fasta'	# Assembly output fasta file

## Merging reference & dark matter DNA
my $hg = 'GRCh38_full_analysis_set_plus_decoy_hla.fa'


## Add dark matter DNA target keyword
my $seqout = Bio::SeqIO->new( -file => '>NRS.fa', -format =>'fasta' );	## output file
my $seqio  = Bio::SeqIO->new( -file => $contigs_in);	## input file
while ( my $seq = $seqio->next_seq ) {
    my $id = $seq->display_id;
    $id =~ s/^(>)/NRS_/;	# Add target keyword
    print $id;	
#    $id =~ s/\_length.*//g;
    $seq->display_id($id);
    $seq->description('');
    $seqout->write_seq($seq);
}

## Merge reference & dark DNA
system( 'cat $hg NRS.fa > GRCh38_nrs.fa' );

## Create bowtie2 index
system( 'bowtie2-build', 'GRCh38_nrs.fa', 'GRCh38_nrs' );
