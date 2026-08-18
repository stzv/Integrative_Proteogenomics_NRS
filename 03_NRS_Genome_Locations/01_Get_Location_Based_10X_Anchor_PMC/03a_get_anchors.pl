#!/usr/bin/perl -w
use strict;

my %need = ();
die "File doesn't exist" unless open F, 'SABE_1172_UNHESMSV_genotyping/SABE1172_UNHESMSV_NRS_dark_freeze_final.fa';
while ( <F> ) {
    $need{$1} = $1.$2 if m/^\>(k141\_\d+)(\_\d+\_\d+)/;
}
close F;
warn "Total NRS: ", scalar keys %need, "\n"; 

my %data = ();
foreach my $file ( <brazil/20210114_genotyping/*_anchors.txt.gz> ) {
    open F, 'gunzip -c '.$file.' |';
    while ( <F> ) {
        my ( $read, $flag, $ctg1, $pos1, $mapq, $cigar, $ctg2, $pos2 ) = split /\t/;
        die 'low mapq: ',$_ if $mapq<20;
        next unless exists( $need{$ctg1} );
        my $strand = '+';
        if ( $ctg1 =~ m/^k141/ and $ctg2 !~ m/^k141/ ) { # R1 - dark contig, R2 - chromosome 
            if ( $flag & 32 ) { # Chromosomal read is minus strand
                $strand = '-' if $flag & 16; # If R1 and R2 both on minus, we will need revcomp contig upon insertion to chr
                push @{$data{$ctg1}}, join( ':', $ctg2, $pos2, '<', $strand );
                
            }
            else { # Chromosomal read is plus strand
                # If a chromosomal read on plus strand then find the genomic coordinates of the end of the read 
                my $pos3 = $pos2 + 149; # Need to find the end of the chromosomal read, dark contig won't start before that :)
                $strand = '-' if not($flag & 16);
                push @{$data{$ctg1}}, join( ':', $ctg2, $pos3, '>', $strand );
            }
        }
    }
    close F;
    warn $file, ' ', scalar keys %data, "\n";
}

open F, '>', '03_SABE_1172_UNHESMSV_anchors.txt';
foreach my $ctg ( sort keys %data ) {
    print F join( "\t", $ctg, sort @{$data{$ctg}} ), "\n";
}
close F;
