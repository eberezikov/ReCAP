#!/usr/bin/perl -w

use strict;

open F,"split.contigs.corrected_orientation.fa";
open OUT,">split.contigs.fa";

local $/ = "\n>";

my $n = 0;
while (my $seq = <F>) {
  $n++;
  $seq=~s/>//g;
  $seq=~s/(.*)\n//;
  my $info = $1;
  $seq=~s/\n//g;
  print OUT ">contig1-$n ",length($seq)," $info\n$seq\n";
}


