#!/usr/bin/perl -w

use strict;

my ($sam1,$sam2,$fq1,$fq2,$parts_fastq,$ctgparts,$K) = @ARGV;

my ($tot,$missed,$mapped_same,$mapped_diff,$mapped_one,$skipped_one,$skipped_both) = (0,0,0,0,0,0,0);
my $n = 0;
my %CTG;
my %HITS1;
open(SAM1, $sam1);
while (my @l = split(/\t/,<SAM1>)) {
   $l[0]=~s/\/.*//;
   push @{$HITS1{$l[0]}},$l[2] unless $l[2] eq '*';
}
close SAM1;

my %HITS2;
open(SAM2, $sam2);
while (my @l = split(/\t/,<SAM2>)) {
   $l[0]=~s/\/.*//;
   push @{$HITS2{$l[0]}},$l[2] unless $l[2] eq '*';
}
close SAM2;

open(FQ1, $fq1);
open(FQ2, $fq2);
while ( my $seq1 = <FQ1>.<FQ1>.<FQ1>.<FQ1>) {
        my $seq2 = <FQ2>.<FQ2>.<FQ2>.<FQ2>;
        my ($id) = $seq1=~/^\@(\S+)\//;
        $tot++;
        unless (defined $HITS1{$id} && defined $HITS2{$id}) {
            $missed++;
            next;
        }
        my %THISCTG;
        my %THISCTG1;
        my $this_skipped = 0;
        my $this_mapped = 0;
        if (defined $HITS1{$id}) {
           unless (scalar @{$HITS1{$id}} == $K) { # skip if max hits reported, since most likely not all are reported
            $this_mapped++;
            foreach my $c (@{$HITS1{$id}}) {
              $THISCTG1{$c}++;
              $THISCTG{$c}++;
            }
           }
           else {
             $this_skipped++;
           } 
         }
         my %THISCTG2;
         if (defined $HITS2{$id}) {
           unless (scalar @{$HITS2{$id}} == $K) { # skip if max hits reported, since most likely not all are reported
            $this_mapped++;
            foreach my $c (@{$HITS2{$id}}) {
              $THISCTG2{$c}++;
              $THISCTG{$c}++;
            }
           }
           else {
             $this_skipped++;
           }
         }
         my $this_same = 0;
         if ($this_skipped == 2) {
            $skipped_both++;
         }
         elsif ($this_skipped == 1) {
           $skipped_one++;
         }
         elsif ($this_mapped == 1) {
           $mapped_one++;
         }
         else {
           foreach my $c (keys %THISCTG1) {
             $this_same++ if defined $THISCTG2{$c};
           }
         }
         if ($this_same > 0) {
           $mapped_same++;
         }
         else {
           $mapped_diff++;
         }
         foreach my $ctg (keys %THISCTG) {
           my ($ctgn) = $ctg=~/\-(\d+)/;
           my $part = int($ctgn/$ctgparts);
           push @{$CTG{$part}}, $ctg."\n".$seq1.$seq2;
         }
}

open OUT,">$parts_fastq";
foreach my $p (sort {$a <=> $b} keys %CTG) {
   print OUT "part$p\n";
   print OUT @{$CTG{$p}};
   print OUT "\n";
}
close OUT;

#     open STAT,">stat.03/$lib.stat.03";
#     print STAT "total pairs:\t$tot\n\n";
#     print STAT "mapped both in a pair to the same contig:\t$mapped_same (",sprintf("%.2f",$mapped_same/$tot*100),"%)\n";
#     print STAT "mapped both in a pair to different contigs:\t$mapped_diff (",sprintf("%.2f",$mapped_diff/$tot*100),"%)\n";
#     print STAT "mapped only one in a pair:\t$mapped_one (",sprintf("%.2f",$mapped_one/$tot*100),"%)\n";
#     print STAT "mapped one in a pair and one is skipped due to multihit:\t$skipped_one (",sprintf("%.2f",$skipped_one/$tot*100),"%)\n";
#     print STAT "skipped both in a pair due to multihit:\t$skipped_one (",sprintf("%.2f",$skipped_both/$tot*100),"%)\n";
#     print STAT "both in a pair are not mapped:\t$missed (",sprintf("%.2f",$missed/$tot*100),"%)\n";
#     my $checksum = $mapped_same+$mapped_diff+$mapped_one+$skipped_one+$skipped_both+$missed;
#     print STAT "\nChecksum: $checksum ",$checksum/$tot,"\n";
#     close STAT;
#     warn "Done with lib $lib\n";


