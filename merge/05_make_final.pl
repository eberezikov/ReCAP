#!/usr/bin/perl -w

use strict;

my $cycle;

if (-e 'MERGING') {
  if (-e 'MERGING/CYCLE_TO_RUN') {
     $cycle = `cat MERGING/CYCLE_TO_RUN`;
     $cycle=~s/\n//;
  }
  else {
    die "No cycle info found\n\n";
  }
}
else {
  die "No MERGING folde exists\n";
}

if ($cycle=~/FINAL/) {
  $cycle=~s/FINAL//;
  $cycle=~s/\n//;
}
else { 
  die "Not a final cycle\n";
}

my $prev = $cycle - 1;

mkdir 'FINAL';

`cp MERGING/c$prev/ctgs_out.nr.fa FINAL/final_all.fa`;
`cp MERGING/c$prev/ctgs_out.nr.fa.qual FINAL/final_all.fa.qual`;

chdir 'FINAL';

warn "Making bowtie index\n";
`/usr/local/bowtie/2.2.4/bowtie2-build final_all.fa final_all.fa 2>log.00 1>&2`;

