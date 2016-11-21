#!/usr/bin/perl -w

use strict;


my @ranks = ([90,1],[95,2],[99,3],[100,4]);


my @ctgs;
my %RANKED;
my %RANK_STAT;

my @header;

foreach my $r (@ranks) {
   my ($cutoff,$rank) = @$r;
   warn "reading $rank $cutoff\n";
   open(F,"FINAL/contigs.cutoff$cutoff.list.classified") || die "No file found!\n";
   @header = split(/\s/,<F>);
   while (my ($id,$len,$type,@libs) = split(/\s/,<F>)) {
      next if defined $RANKED{$id};
      push @ctgs,[$id,$len,$rank,$type,@libs];
      $RANKED{$id}=$rank;
      $RANK_STAT{$rank}++;
   }
   close F;
}



open OUT,">FINAL/contigs.ranked";
print OUT join("\t",'id','len','rank','type',@header[3..$#header]),"\n";
foreach my $c (@ctgs) {
   print OUT join("\t",@$c),"\n";
}
close OUT;

open OUT,">FINAL/contigs.ranked.stats";
foreach my $rank (sort {$a <=> $b} keys %RANK_STAT) {
   print OUT $rank,"\t",$RANK_STAT{$rank},"\n";
}












