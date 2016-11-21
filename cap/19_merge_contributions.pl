#!/usr/bin/perl -w

use strict;

#my @contrib_files = @ARGV;
my @contrib_files = <FINAL/map/*/contributions.txt>;

my $MINREADS = 1;
my $MIN_FRACTION_FOR_CLUSTERING = 0.5;

my %TR;
my %TRREADS;

foreach my $f (@contrib_files) {
 open(F,$f) || die "Can't open $f\n";
 warn "Reading $f\n";
 while (my ($id,@trs) = split(/\s/,<F>)) {
   foreach my $t (@trs) {
      my ($id1,$r) = split(/\:/,$t);
      if ($r >= $MINREADS) {
         $TR{$id}{$id1}+=$r;
         $TRREADS{$id}+=$r;
      }
   }
 }
}
close F;


open OUT,">FINAL/contributions_merged.txt";
foreach my $id (sort {$TRREADS{$b} <=> $TRREADS{$a}} keys %TRREADS) {
    my $totreads = 0;
    my $selfreads = 0;
    my @contributors = sort {$TR{$id}{$b} <=> $TR{$id}{$a}} keys %{$TR{$id}};
    foreach my $id1 (@contributors) {
       $totreads+=$TR{$id}{$id1};
       if ($id1 eq $id) {
         $selfreads = $TR{$id}{$id1};
       }
    }
    print OUT join("\t",$id,$totreads,$selfreads,scalar @contributors);
    foreach my $c (@contributors) {
      print OUT "\t$c:$TR{$id}{$c}";
    }
    print OUT "\n";
}

