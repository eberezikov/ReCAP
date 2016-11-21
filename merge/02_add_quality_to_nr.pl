#!/usr/bin/perl -w

use strict;

my %QUAL;
open(F,"assembled_ctgs.00.fa.qual");
while (my $h = <F>) {
  my $qual = <F>;
  my ($id) = $h=~/>(\S+)/;
  $QUAL{$1}=$h.$qual;;
}
close F;

open(F,"assembled_ctgs.00.nr99.fa");
open(OUT,">assembled_ctgs.00.nr99.fa.qual");
while (my $h = <F>) {
  my $seq = <F>;
  my ($id) = $h=~/>(\S+)/;
  print OUT $QUAL{$id};
}
