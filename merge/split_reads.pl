#!/usr/bin/perl -w

use strict;
use Parallel::ForkManager;

my ($indir, $outdir, $seqs_per_part) = @ARGV;

$seqs_per_part = 100000 unless defined $seqs_per_part;

my $pm = Parallel::ForkManager->new(8);

unless (-e $indir) {
  die "No $indir found\n";
}

if (-e $outdir) {
#  die "$outdir already exists!\n";
}
else {
  mkdir $outdir;
}


while (my $f1 = <$indir/*/*.R1.fastq>) {
 $pm->start and next;
  my $f2 = $f1; $f2=~s/\.R1\./\.R2\./;
  my ($lib) = $f1=~/$indir\/(\S+)\//;
  my $libdir = "$outdir/$lib";
  if (-e $libdir) { die "\n\n!!! $libdir already existst\n\n"; }
  warn "$f1\n$f2\n$lib\n\n";
#  open(R1,"gzip -dc $f1 |");
#  open(R2,"gzip -dc $f2 |");
  open(R1,"$f1");
  open(R2,"$f2");
  my $part = 1;
  my $curr_reads = 0;
  mkdir "$libdir";
  open OUTR1, ">$libdir/$lib.p$part.R1.fq";
  open OUTR2, ">$libdir/$lib.p$part.R2.fq";
  warn "$lib $part\n";
  my $n = 0;
  while (my $h1 = <R1>) {
    my $seq1 = <R1>;
    my $hq1 = <R1>;
    my $qual1 = <R1>;
    my $h2 = <R2>;
    my $seq2 = <R2>;
    my $hq2 = <R2>;
    my $qual2 = <R2>;
    $curr_reads++;
    $n++;
    if ($curr_reads > $seqs_per_part) {
      close OUTR1;
      close OUTR2;
      $part++;
      open OUTR1, ">$libdir/$lib.p$part.R1.fq";
      open OUTR2, ">$libdir/$lib.p$part.R2.fq";
      $curr_reads = 1;
      warn "$lib $part\n";
    }
    print OUTR1 $h1,$seq1,$hq1,$qual1;
    print OUTR2 $h2,$seq2,$hq2,$qual2;
  }
  close OUTR1;
  close OUTR2;
 $pm->finish;
}

$pm->wait_all_children();