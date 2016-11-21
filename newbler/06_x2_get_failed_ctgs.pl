#!/usr/bin/perl -w

use strict;
use Config::General qw(ParseConfig);
my %CONF = ParseConfig("CONFIG");


my $ctgparts = $CONF{'CTGPARTS'};


if (-e "ctgs_rerun") {
   die "ctgs_rerun exists!\n";
}
else {
  mkdir "ctgs_rerun";
}

my %CTG;
my %PART;

open F,"failed.06";
while (my $l = <F>) {
   my ($reads,$id) = $l=~/with (\d+) reads for (contig\S+)\!/;
#   next if $reads < 100;
   $CTG{$id}++;
   my ($part) = $id=~/-(\d+)/;
   $part=int($part/$ctgparts);
   $PART{$part}++;
}
close F;

my @to_parse;
while (my $f = <ctgs_large/*.ctg>) {
  my ($id) = $f=~/\/(contig\S+)\.ctg/;
  push @to_parse, $f if defined $CTG{$id};
}

while (my $f = <ctgs_small/*.ctg>) {
  my ($p) = $f=~/part(\d+)/;
  push @to_parse, $f if defined $PART{$p};
}


my $n = 0;
my $tot = scalar keys %CTG;
foreach my $f (@to_parse) {
  open(F,$f);
  warn "$f\n";
  local $/ = "\n\n";
  while (my $block=<F>) {
     my ($id) = $block=~/^(\S+)/;
     next unless defined $CTG{$id};
     $n++;
     warn "$n/$tot. $id. $f\n";
     open(OUT,">ctgs_rerun/$id.ctg");
     print OUT $block;
     close OUT;
  }
  close F;
}
