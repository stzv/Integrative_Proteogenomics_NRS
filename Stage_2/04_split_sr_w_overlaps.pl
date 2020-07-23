#!/usr/bin/perl -w
use strict;
## Creats 3 sam files with mapped, unmapped and partially mapped as well as map.txt file with partially mapped (0/1)
my $min_disc = 0.05;	# 95% Identity match
my $min_unmap = 200;	# Minimal length of unmapped NRS section

my $minimap_alignment = 'NRS_vs_GRCh38.sam.gz';
my $contig_length_file = 'contig_lens.txt.gz';

my %mapped = ();
open F, 'gunzip -c '.$minimap_alignment.' |';
open F1, '| gzip -c >04_unmapped.sam.gz';
open F2, '| gzip -c >04_partial.sam.gz';
open F3, '| gzip -c >04_mapped.sam.gz';
while ( <F> ) {
    next if m/^\@/;
    my ( $ctg, $flag, $chr, $pos, $mapq, $cigar ) = split /\t/;
    next if $flag & 2048; 
    my $div = 1;
    if ( m/de\:f\:([\d\.]+)/ ) { 
        $div = $1;
    }
    if ( $div > $min_disc ) { 
        print F1;
        $mapped{$ctg} = 1;
    }
    elsif ( $cigar =~ m/(\d{3,})[SDI]/ and $1 >= $min_unmap ) { 
        print F2;
        $mapped{$ctg} = 2;
    }
    else { 
        print F3;
        $mapped{$ctg} = 3;
    }
}
close F;
my %maps = ();
my %coord = (); 
warn "Reading Contig lengths\n"; 
open F, 'gunzip -c '.$contig_length_file.' |'; 
while ( <F> ) {
    chomp;
    next unless m/^\@SQ\tSN\:(NRS\d+)\tLN\:(\d+)/; 
    my ( $ctg, $len ) = ( $1, $2 ); 
    $maps{$ctg} = '0'x$len if $mapped{$ctg} == 2; 
}
close F;
open F, 'gunzip -c '.$minimap_alignment.' |';
while ( <F> ) {
    next if m/^\@/;
    my ( $ctg, $flag, $chr, $pos, $mapq, $cigar ) = split /\t/;
    #next if $flag & 2048;
    next unless exists($maps{$ctg}); 
    $coord{$ctg} = join( ':', $flag, $chr, $pos, $mapq, $cigar ) unless $flag & 2048; 
    my $div = 1;
    if ( m/de\:f\:([\d\.]+)/ ) { 
        $div = $1;
    }
    next if $div > $min_disc; 
    my $cur = 0;
    my $map = $maps{$ctg};
    $map = reverse($map) if $flag & 16; 
    while ( $cigar =~ m/(\d+)(\D)/g ) {  
        my ( $len, $type ) = ( $1, $2 );
        if ( $type eq 'M' or $type eq 'I' ) {
            substr($map, $cur, $len, '1'x$len);
        }
        $cur += $len unless $type eq 'D';
    }
    $map = reverse($map) if $flag & 16;
    $maps{$ctg} = $map;
}
close F;

my ($i,$j) = (0,0);
warn "Saving data\n";
open F, '| gzip -c >04_map.txt.gz';
foreach my $ctg( sort keys %maps ) {
    if ( $maps{$ctg} =~ m/^(0{$min_unmap,})/ or $maps{$ctg} =~ m/(0{$min_unmap,})$/ ) {
        print F $ctg, "\t", $coord{$ctg}, "\t", $maps{$ctg}, "\n";
        $i++;
        $j+=length($1);
    } 
}
close F;
warn "Total partial: $i bp:$j out of ",scalar(keys %maps),"\n";
