#!/usr/bin/perl -w

use strict;

# 0       1       2                3                4              5                6              7               8               9               10     11      12       13             14               15               16              17         18               19               20               21              22
#ID      len     DV1polyflex     DV1riboflex     Low_15h_3       Low_15h_4       Low_1h_1        Low_1h_2        Low_60h_5       Low_60h_6       mlsseq2 mlsseq3 mlsseq4 RAMPAGE1        RAMPAGE2        Starved_15h_3   Starved_15h_4   Starved_1h_1 Starved_1h_2    Starved_60h_5   Starved_60h_6   Stijn20worms    Stijn2Ug        Totreads

my %LIBTYPE = ('cds' => ['DV1polyflex','Low_15h_3','Low_15h_4','Low_1h_1','Low_1h_2','Low_60h_5','Low_60h_6','Starved_15h_3','Starved_15h_4','Starved_1h_1','Starved_1h_2','Starved_60h_5','Starved_60h_6','Stijn20worms','Stijn2Ug'],
               'utr3' => ['mlsseq2','mlsseq3','mlsseq4'],
               'utr5' => ['RAMPAGE1','RAMPAGE2'],
               'nc' => ['DV1riboflex'],
              );

my $min_reads = 3;

foreach my $t (keys %LIBTYPE) {
   open OUT,">FINAL/$t.counts";
   my %CTG;
   foreach my $f (@{$LIBTYPE{$t}}) {
     $f="FINAL/map/$f/read_counts.txt";
     warn "readig for $t from $f\n";
     open F,"$f";
     while (my ($ctg,$reads) = split(/\s/,<F>)) {
       $CTG{$ctg}+=$reads;
     }
     close F;
   }
   foreach my $ctg (sort keys %CTG) {
      print OUT $ctg,"\t",$CTG{$ctg},"\n";
   }
   close OUT;
}
