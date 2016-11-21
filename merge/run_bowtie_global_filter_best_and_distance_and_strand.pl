#!/usr/bin/perl -w

use strict;

my ($genome,$fq,$sam,$err,$P,$K,$max_distance,$strand) = @ARGV;

my $prev_id = "";
my $previous_score = 0;
$max_distance = 5 unless defined $max_distance;
$strand = 'F' unless defined $strand;

open (BWT,"/usr/local/bowtie/2.2.4/bowtie2 -p $P -k $K -x $genome -q -U $fq --end-to-end  --reorder --sam-nohead --sam-nosq 2>$err |");

open OUT,">$sam";
#open OUT1,">$sam.orig";

while (my $l = <BWT>) {
#  print OUT1 $l;
  my ($id,$flag) = $l=~/^(\S+)\t(\d+)/;
  next if ($flag == 4); #unpammed
  if ($flag == 0 || $flag == 256) { #plus-mapped
    next if ($strand ne 'F');
  } 
  elsif ($flag == 16 || $flag == 272) {
    next if ($strand ne 'R');
  }
  if ($l=~/AS\:i\:(-?\d+).*\tNM\:i\:(\d+)/) {
       my ($score,$distance) = ($1,$2);
       next if ($distance > $max_distance);
       next if (($id eq $prev_id) && ($score < $previous_score));
       $prev_id = $id;
       $previous_score = $score;
       print OUT $l;
  }
}

