#!/usr/bin/perl -w
use strict;

my $min_cov = 7.5;
my $max_cov = 100;

my %lens = ();
open F, 'gunzip -c contig_lens.txt.gz |';
while ( <F> ) {
    next unless m/^\@SQ\tSN\:(.*)\tLN\:(\d+)/;
    $lens{$1} = $2; 
}
close F;
warn "Total:",scalar keys %lens, "\n";

my $total = 0;
my %non38 = ();
open F, 'gunzip -c 04_unmapped.sam.gz |';
while ( <F> ) {
    my ( $ctg ) = split /\t/;
    $non38{$ctg} = 1; 
    $total += $lens{$ctg};    
}
close F;
warn "Total: $total bases in ",scalar keys %non38, " contigs\n";

my %lens2 = (); 
open F, 'gunzip -c 04_map.txt.gz |'; 
while ( <F> ) {
    chomp;
    my ( $ctg, $mapping, $bitstr ) = split /\t/;
    $non38{$ctg} = 2;    
    $total += $lens{$ctg}; 
    if ( $ctg =~ m/^1+.*?(0{200,})$/ ) { 
        $lens2{$ctg} =length($1);   
    }
    elsif ( $ctg =~ m/^(0{200,}).*?1+$/ ) { 
        $lens2{$ctg} =length($1);   
    }
    else {
        $lens2{$ctg} = $lens{$ctg};   
    }
}
close F;
warn "Total: $total bases in ",scalar keys %non38, " contigs\n";
my %table = ();  
my @samples = ();
my $total_len = 0;
foreach my $file ( <bowtie_vs_extended/*_dark_cnt.txt.gz> ) {
    my $sample = $file;
    $sample =~ s/.*\///g;
    $sample =~ s/\_unmapped\_dark_cnt\.txt\.gz$//;
    warn "Reading $sample\n";
    push @samples, $sample;
    open F, 'gunzip -c '.$file.' |';
    my $new = 0;
    while ( <F> ) {
        chomp;
        my ( $ctg, $count ) = split /\t/;
        next unless exists($non38{$ctg});
        $table{$ctg}{$sample} = $count;
#        $new++ if $count >= 10 and !exists($table{$ctg}{'total'});
#        $table{$ctg}{'total'}+=$count if $count >=10;
        if ( $count * 150 >= $min_cov * $lens{$ctg} and $count * 150 <= $max_cov * $lens{$ctg} and !exists($table{$ctg}{'total'}) ) { 
            $total_len += exists( $lens2{$ctg} ) ? $lens2{$ctg} : $lens{$ctg};
            $new++; ##just a counter
        }
        $table{$ctg}{'total'}+=$count if $count * 150 >= $min_cov * $lens{$ctg} and $count * 150 <= $max_cov * $lens{$ctg}; 
    }
    close F;
    warn $file, ' ', scalar( keys %table), " new: $new, total bases: $total_len\n";
}
my ( $total1, $total2, $bp1, $bp2 ) = ( 0, 0, 0, 0 ); 
open F, '>', '05_freq_results.txt';
print F join( "\t", 'Contig', @samples ), "\n";
foreach my $ctg ( sort keys %table ) {
#    next if $table{$ctg}{'total'} < 10;
    next unless exists( $table{$ctg}{'total'} );
    $bp1 += $lens{$ctg} if $non38{$ctg} == 1;
    $total1++ if $non38{$ctg} == 1;
    $bp2 += $lens2{$ctg} if $non38{$ctg} == 2;
    $total2++ if $non38{$ctg} == 2;
    print F $ctg;
    foreach my $sample ( @samples ) {
        my $cnt = exists( $table{$ctg}{$sample} ) ? $table{$ctg}{$sample} : 0; #what is written in the file?
        print F "\t", $cnt;
    }
    print F "\n";
}
close F;
warn "Fully unmapped: $total1, $bp1 bases\n";
warn "Partially unmapped: $total2, $bp2 bases\n";