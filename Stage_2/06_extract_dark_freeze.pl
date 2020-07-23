#!/usr/bin/perl -w
use strict;
use lib '/usr/local/ensembl72/bioperl-live';
use Bio::SeqIO;

my $infile = 'NRS.fa';
my $outfile = 'NRS_dark_freeze.fa';

my %need = (); 
open F, '05_freq_results.txt';
<F>; 
while ( <F> ) {
    my ( $ctg ) = split /\t/;
    $need{$ctg} = 1; 
}
close F;
warn "Need ", scalar keys %need, "\n";

my %partial = ();
my %placing = ();
open F, 'gunzip -c 04_map.txt.gz |';
while ( <F> ) {
    chomp;
    my ( $ctg, $info, $bitstr ) = split /\t/;
    next unless exists($need{$ctg}); 
    if ( $bitstr =~ m/^0{200,}/ and $bitstr =~ m/0{200,}$/ and $bitstr !~ m/^0+$/ ) { 
        $partial{$ctg} = '1_'.length($bitstr); 
    }  
    elsif ( $bitstr =~ m/^(0{200,})/ ) { 
        $partial{$ctg} = '1_'.length($1);
        $placing{$ctg} = $info;
    }
    elsif ( $bitstr =~ m/(0{200,})$/ ) {
        $partial{$ctg} = length($`).'_'.length($bitstr); 
        $placing{$ctg} = $info;
    }
    else {
        die '11 '.$ctg; 
    }
}
close F;
warn scalar keys %partial, " partial processed\n";

my %ctg2id = ();
my ( $max_len, $tot_len ) = ( 0, 0 );
my $seqio = Bio::SeqIO->new( -file => $infile, -format => 'fasta' );
open F, '>', $outfile;
while ( my $seq = $seqio->next_seq ) {
    my $ctg = $seq->display_id;
    next unless exists( $need{$ctg} );
    my $seqstr = $seq->seq;
    my ( $s, $e ) = ( 1, length($seqstr) );
    $max_len = $e if $max_len < $e;
    if ( exists( $partial{$ctg} )) {
        ( $s, $e )= split /\_/, $partial{$ctg};
        $seqstr = substr( $seqstr, $s-1, $e-$s+1 ); 
    }
    $tot_len += length($seqstr);
    $seqstr =~ s/(\w{60})/$1\n/g;
    $seqstr =~ s/\n$//;
    $ctg2id{$ctg}=$ctg.'_'.$s.'_'.$e;
    if ( exists($placing{$ctg}) ) { 
        print F '>', $ctg, '_',$s, '_',$e, ' ', $placing{$ctg}, "\n", $seqstr, "\n";
    }
    else {
        print F '>', $ctg, '_',$s, '_',$e, "\n", $seqstr, "\n";
    }
}
close F;
warn "Max:",$max_len, "\n";
warn "Tot:",$tot_len, "\n";

open F, '05_freq_results.txt';
open F1, '>', 'NRS_dark_freeze_coverage.txt';
my $header = <F>;
print F1 $header;
while( <F> ) {
    chomp;
    my @arr = split /\t/;
    next unless exists( $ctg2id{$arr[0]} );
    print F1 join( "\t", $ctg2id{$arr[0]}, @arr[1..$#arr] ), "\n";
}
close F;
close F1;
