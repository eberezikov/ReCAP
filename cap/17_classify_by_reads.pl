#!/usr/bin/perl -w

use strict;

# 0       1       2                3                4              5                6              7               8               9               10     11      12       13             14               15               16              17         18               19               20               21              22
#ID      len     DV1polyflex     DV1riboflex     Low_15h_3       Low_15h_4       Low_1h_1        Low_1h_2        Low_60h_5       Low_60h_6       mlsseq2 mlsseq3 mlsseq4 RAMPAGE1        RAMPAGE2        Starved_15h_3   Starved_15h_4   Starved_1h_1 Starved_1h_2    Starved_60h_5   Starved_60h_6   Stijn20worms    Stijn2Ug        Totreads

my @cds = (2,4,5,6,7,8,9,15,16,17,18,19,20,21);
my @utr3 = (10,11,12);
my @utr5 = (14,15);
my @nc = (3);

my $min_reads = 3;

while (my $f = <FINAL/contigs.cutoff*.list>) {
   warn "$f\n";
   open F,$f;
   open OUT,">$f.classified\n";
   <F>;
   my %TYPES;
   print OUT join("\t",'id','len','type','utr5','cds','utr3','nc'),"\n";
   while (my @l = split(/\s/,<F>)) {
      my $id = $l[0];
      my $len = $l[1];
      my $tot_cds = 0;
      for my $i (@cds) {
         $tot_cds+=$l[$i];
      }
      my $tot_utr3=0;
      for my $i (@utr3) {
         $tot_utr3+=$l[$i];
      }
      my $tot_utr5=0;
      for my $i (@utr5) {
         $tot_utr5+=$l[$i];
      }
      my $tot_nc=0;
      for my $i (@nc) {
         $tot_nc+=$l[$i];
      }

      my @type;
      if ($tot_utr5 >= $min_reads) {
         push @type, '5UTR';
      }
      if ($tot_cds >= $min_reads) {
         push @type, 'CDS';
      }
      if ($tot_nc >= $min_reads && $tot_nc > 2*$tot_cds) {
         push @type, 'NC';
      }
      if ($tot_utr3 >= $min_reads) {
         push @type, '3UTR';
      }
      
      push @type, 'LOW' unless @type;

      my $typestr = join('_',@type);
      $TYPES{$typestr}++;

      print OUT join("\t",$id,$len,$typestr,$tot_utr5,$tot_cds,$tot_utr3,$tot_nc),"\n";
   }
   close F;
   close OUT;
   open OUT,">$f.classified.stats";
   foreach my $type (sort {$TYPES{$b} <=> $TYPES{$a}} keys %TYPES) {
      print OUT $type,"\t",$TYPES{$type},"\n";
   }
   close F;
   close OUT;

}


