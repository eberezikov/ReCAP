#!/usr/bin/perl -w


use strict;

my ($newcycle,$ctgsin,$ctgsout) = @ARGV;

open(F,$ctgsin);
open(Q,"$ctgsin.qual");
open(OUT,">$ctgsout");
open(OUTQUAL,">$ctgsout.qual");

while (my $h = <F>) {
  my $seq = <F>;
  my $hq = <Q>;
  my $qual = <Q>;
  my ($n) = $h=~/>cap(\d+)/;
  if ($n == $newcycle) {
    print OUT $h,$seq;
    print OUTQUAL $hq,$qual;
  }
}

close F;
close Q;
close OUT;
close OUTQUAL;

