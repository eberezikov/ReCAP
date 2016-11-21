#!/usr/bin/perl -w

use strict;

my ($in_pairs,$out_pairs,$max_linked_ctgs) = @ARGV;

$max_linked_ctgs = 5 unless defined $max_linked_ctgs;

my %CTG;
my @pairs;

open(F,$in_pairs);
while (my @c = split(/\s/,<F>)) {
  $CTG{$c[0]}++;
  $CTG{$c[1]}++;
  push @pairs,[@c];
  
}
close F;

my %USED_CTGS;

open OUT,">$out_pairs";
foreach my $p (@pairs) {
   next if $CTG{$p->[0]} > $max_linked_ctgs;
   next if $CTG{$p->[1]} > $max_linked_ctgs;
   next if defined $USED_CTGS{$p->[0]};
   next if defined $USED_CTGS{$p->[1]};
   print OUT join("\t",@$p),"\n";
   $USED_CTGS{$p->[0]}++;
   $USED_CTGS{$p->[1]}++;
}