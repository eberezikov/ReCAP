#!/usr/bin/perl -w

use strict;

my @cutoffs = (90,95,97,99,100);

my $totreads;
my @ctg;
my %TOT_CTG;
open F,"FINAL/clusters.txt";
warn "reading clusters\n";
while (my ($cl,$reads,$n,@members) = split(/\s/,<F>)) {
     foreach my $m (@members) {
        my ($id,$r) = split(/\:/,$m);
        push @ctg,[$cl,$id,$r];
        $totreads+=$r;
        $TOT_CTG{$id}++;
     }
}
close F;

my $totnonzeroctg = scalar keys %TOT_CTG;

warn "sorting clusters\n";
@ctg = sort {$b->[2] <=> $a->[2] || $a->[1] cmp $b->[1]} @ctg;

foreach my $cutoff (@cutoffs) {
      warn "Generating cut-off $cutoff\n";
      my $cutoff_reads = int($totreads*$cutoff/100)+1;
      my %CL;
      my %CLREADS;
      my $cumulative=0;
      my %USED_CTGS ;
      foreach my $c (@ctg) {
         $cumulative+=$c->[2];
         push @{$CL{$c->[0]}},$c->[1].':'.$c->[2];
         $CLREADS{$c->[0]}+=$c->[2];
         $USED_CTGS{$c->[1]}++;
         last if $cumulative >= $cutoff_reads;
      }
      my $totctg = scalar keys %USED_CTGS;
      my $totclusters = scalar keys %CL;
      open OUT,">FINAL/clusters.cutoff$cutoff.txt";
      print OUT "Total reads: $totreads\n";
      print OUT "Total nonzero contigs: $totnonzeroctg\n";
      print OUT "Total used contigs: $totctg\n";
      print OUT "Total generated clusters: $totclusters\n";
      print OUT "Reads for $cutoff\% cutoff: $cutoff_reads\n\n";
      
      foreach my $c (sort {$CLREADS{$b} <=> $CLREADS{$a}} keys %CLREADS) {
        print OUT join("\t",$c,$CLREADS{$c},scalar @{$CL{$c}},@{$CL{$c}}),"\n";
      }

      close OUT;
}
