#!/usr/bin/perl -w
use strict;

## ENVIRONMENT VARIABLES ##
my $min_mapq = 20;								# minimal mapping quality
my $nrs_location = '/mnt/fedot21/sabe/*.fq.gz'	# location folder of the FastQ sample files
my $threads = 48								# Number of threads

##Opening the initial fastq file and deviding it into two 1.fq and 2.fq files (pair-end)
foreach my $fq ( $nrs_location ) {
    my $sample = $fq;
    $sample =~ s/.*\///g;
    $sample =~ s/\_.*$//g;
    next if -e $sample.'_dark_cnt.txt.gz';
    warn "Preparing $sample\n";
    open F, 'gunzip -c '.$fq.' |';
    open F1, '>', '1.fq';
    open F2, '>', '2.fq';
    while ( my $id1 = <F> ) {
        my $seq1 = <F>;
        my $plus1 = <F>;
        my $qual1 = <F>;
        my $id2 = <F>;
        my $seq2 = <F>;
        my $plus2 = <F>;
        my $qual2 = <F>;
        die "ID conflict $id1 $id2" unless $id1 eq $id2;
        print F1 $id1,$seq1,$plus1,$qual1;
        print F2 $id2,$seq2,$plus2,$qual2;
    } 
    close F1;
    close F2;
    close F;
    warn "Mapping $sample\n";
    my %counts = ();
    open F, join( ' ', 'bowtie2',
                              '-x', 'GRCh38_nrs',
                              '-1', '1.fq',
                              '-2', '2.fq',
                              '-p', $threads,
                              '|' );
    open F1, '| gzip -c >bowtie_vs_extended/'.$sample.'_anchors.txt.gz';
    while ( <F> ) {
        next if m/^\@/;
        my ( $read, $flag, $chr, $pos, $mapq, $cigar, $chr2, $pos2, $tlen ) = split /\t/;
        if ( $chr =~ m/^NRS/ and $mapq >= $min_mapq ) {
            $counts{$chr}++;
            print F1 if $chr2 ne '*' and $chr2 ne '='; 
        }
        $counts{$chr}++ if $chr =~ m/^NRS/ and $mapq >= $min_mapq;
    }
    close F;
    warn "Outputting $sample\n";
    open F, '| gzip -c >bowtie_vs_extended/'.$sample.'_dark_cnt.txt.gz';
    foreach my $ctg (sort keys %counts) {
        print F $ctg, "\t", $counts{$ctg}, "\n";
    }
    close F;
}