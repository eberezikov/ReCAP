#!/usr/bin/perl -w

use strict;

my @cutoffs = (90,95,97,99,100);
my @libs;
while (my $l = <FINAL/map/*>) {
   $l=~s/.*\///;
   push @libs, $l;
}

open(STAT,">FINAL/contigs.cutoff_stats");

foreach my $cutoff (@cutoffs) {
   my %CTG;
   foreach my $lib (@libs) {
      next unless -e "FINAL/map/$lib/cutoff$cutoff.txt";
      open(F,"FINAL/map/$lib/cutoff$cutoff.txt");
      warn "Cutoff $cutoff: processing $lib\n";
      <F>; <F>; <F>; <F>; <F>;
      while (my @l = split(/\s/,<F>)) {
        $CTG{$l[0]}{'len'} = $l[1];
        $CTG{$l[0]}{$lib} = $l[2];
      }
      close F;
   }
   print STAT "Cutoff $cutoff: ",scalar keys %CTG,"  contigs\n";
   warn "Cutoff $cutoff: writing results\n";
   open(OUT,">FINAL/contigs.cutoff$cutoff.list");
   print OUT join("\t",'ID','len',@libs,'Totreads'),"\n";
   my @data;
   foreach my $id (keys %CTG) {
      my @ctg;
      my $sum = 0;
      push @ctg,$id,$CTG{$id}{'len'};
      foreach my $lib (@libs) {
         $CTG{$id}{$lib}=0 unless defined $CTG{$id}{$lib};
         push @ctg, $CTG{$id}{$lib};
         $sum+=$CTG{$id}{$lib};
      }
      push @data,[$sum,@ctg];
   }
   foreach my $d (sort {$b->[0] <=> $a->[0]} @data) {
      my ($sum,@res) = @$d;
      print OUT join("\t",@res,$sum),"\n";
   }
#   last;
}
