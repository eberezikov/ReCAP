#!/usr/bin/perl -w

use strict;

my ($workdir,$cycle) = @ARGV;

chdir("$workdir/MERGING/c$cycle") || die "Can't go to $workdir/MERGING/c$cycle\n";

$cycle++;

my $new_merges = 0;

my %QUAL;
open(F,"ctgs_out.all.fa.qual");
while (my $h = <F>) {
  my $qual = <F>;
  my ($id) = $h=~/>(\S+)/;
  $QUAL{$id}=$h.$qual;;
}
close F;

open(F,"ctgs_out.nr.fa");
open(OUT,">ctgs_out.nr.fa.qual");
while (my $h = <F>) {
  my $seq = <F>;
  my ($id) = $h=~/>(\S+)/;
  print OUT $QUAL{$id};
  if ($id=~/^cap$cycle/) { $new_merges++; }
}
close OUT;

chdir('..');

open OUT,">CYCLE_TO_RUN";
if ($new_merges > 0) {
  print OUT $cycle,"\n";
}
else {
  print OUT "FINAL\n$cycle\n";
}
close OUT;

