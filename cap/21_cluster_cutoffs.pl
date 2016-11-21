#!/usr/bin/perl -w

use strict;

my @cutoffs = (90,95,97,99,100);

my $tot;
my @ctg;
open F,"FINAL/clusters.txt";
while (my @l = split(/\s/,<F>)) {
      push @ctg,[@l];
      $tot+=$l[1];
}
close F;


foreach my $cutoff (@cutoffs) {
      my $cutoff_reads = int($tot*$cutoff/100)+1;
      my $totctg = scalar @ctg;
      my @passed;
      open OUT,">FINAL/clusters.cutoff$cutoff.txt";
      print OUT "Total reads: $tot\n";
      print OUT "Total non-zero contigs: $totctg\n";
      print OUT "Reads for $cutoff\% cutoff: $cutoff_reads\n";

      my $cumulative = 0;
      foreach my $c (@ctg) {
        $cumulative+=$c->[1];
        push @passed, $c;
        last if $cumulative >= $cutoff_reads;
      }

      print OUT "Contigs at $cutoff\% cutoff: ", scalar @passed,"\n\n";

      foreach my $c (@passed) {
        print OUT join("\t",@$c),"\n";
      }

      close OUT;
}
