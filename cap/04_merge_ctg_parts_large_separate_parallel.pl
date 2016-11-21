#!/usr/bin/perl -w

use strict;
use Parallel::ForkManager;

if (-e "ctgs_small") { die "ctgs_small already exists!\n"; }
else { mkdir "ctgs_small" }

if (-e "ctgs_large") { die "ctgs_large already exists!\n"; }
else { mkdir "ctgs_large" }

my $max_reads = 50000;

my %PARTS;

while (my $ctgfile = <LIBS_LS/*/*.ctg>) {
   my ($part) = $ctgfile=~/part(\d+)\.ctg/;
   push @{$PARTS{$part}},$ctgfile;
}

my %CTGSEQ;
for my $f ('split.contigs.fa') {
  warn "Reading $f\n";
  local $/ = "\n>";
  open F, $f;
  while (my $seq = <F>) {
    $seq=~s/>//g;
    $seq=~s/(.*)\n//;
    my $id = $1;
    $id=~s/\s.*//;
    $seq=~s/\n//g;
    $CTGSEQ{$id}=$seq;
  }
  close F;
}

local $/ = "\n";

my $pm = Parallel::ForkManager->new(6);

foreach my $p (sort {$a <=> $b} keys %PARTS) {
  $pm->start and next;
     warn "Merging part $p\n";
     my %CTG;
     foreach my $f (@{$PARTS{$p}}) {
       open F, $f;
       while (my $ctg = <F>) {
         $ctg=~s/\n//;
         my $seqs = <F>.<F>.<F>.<F>.<F>.<F>.<F>.<F>;
         push @{$CTG{$ctg}}, $seqs;
       }
       close F;
     }
     open OUT, ">ctgs_small/part$p.ctg";
     foreach my $ctg (sort keys %CTG) {
       if (scalar @{$CTG{$ctg}} > $max_reads) {
          open OUTLARGE,">ctgs_large/$ctg.ctg";
          print OUTLARGE $ctg,"\n",$CTGSEQ{$ctg},"\n",'0 'x(length($CTGSEQ{$ctg})),"\n",@{$CTG{$ctg}},"\n";
          close OUTLARGE;
       }
       else {
         print OUT $ctg,"\n",$CTGSEQ{$ctg},"\n",'0 'x(length($CTGSEQ{$ctg})),"\n",@{$CTG{$ctg}},"\n";
       }
     } 
     close OUT;
  $pm->finish;
}

$pm->wait_all_children();

