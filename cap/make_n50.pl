#!/usr/bin/perl -w

use strict;

my @ctg;

my $totlen=0;

open F, shift;
local $/ = "\n>";
while (my $seq = <F>) {
  $seq=~s/>//g;
#  next unless $seq=~/^(ctg)/;
  $seq=~s/(.*)\n//;
  my $id = $1;
  $seq=~s/\n//g;
  my $len = length($seq);
#  if ($len < 300) {
#    warn "$id:\n$seq\n";
#  }
  push @ctg, $len;
  $totlen+=$len;
}

@ctg = sort {$b <=> $a} @ctg;

my $n50bases = $totlen/2;

my $curr_len = 0;
my $curr_n = 0;

my $totctg = scalar @ctg;
print "Total contigs: $totctg\n";
print "Total len: $totlen\n";
print "Average len: ",$totlen/$totctg,"\n";
print "Shortest: ",$ctg[$#ctg],"\n";
print "Longest: ",$ctg[0],"\n";

for my $c (@ctg) {
  $curr_len+=$c;
  $curr_n++;
  if ($curr_len >= $n50bases) {
     print "N50: $c ($curr_n scaffolds)\n";
     last;
 }
}
