#!/usr/bin/perl -w

use strict;

my ($libnames,$input_db,$max_distance,$max_freq) = @ARGV;

my %SIZES;
my %ORIENTATION;
open(F,$libnames) || die "No LIBNAMES\n";
while (my @l = split(/\s/,<F>)) {
   $SIZES{$l[1]}{'min'}=$l[2];
   $SIZES{$l[1]}{'max'}=$l[3];
   $ORIENTATION{$l[1]}=$l[4];
}
close F;

my %CTGSIZE;
open(F,$input_db);
while (my $h = <F>) {
  my ($id,$size) = $h=~/>(\S+) (\d+)/;
  $CTGSIZE{$id}=$size;
  <F>;
}
close F;

$max_distance = 3 unless defined $max_distance;

open PAIRS,">pairs";
open LINKS,">links";
my %MAPPED;
my %FREQ;
my %USED_READS;
while (my $f1 = <*.R1.sam>) {
   my %READS1 = %{read_multisam($f1)};
   my $f2=$f1;
      $f2=~s/R1\.sam/R2\.sam/;
   my %READS2 = %{read_multisam($f2)};
   foreach my $id (sort keys %READS1) {
#      warn "checking $id: ", join(" ",@{$READS1{$id}}),"\n";
      if (defined $READS2{$id}) {
#           warn " has a pair: ", join(" ",@{$READS2{$id}}),"\n";
           my ($lib) = $id=~/^(L\d+)/;
           my ($r) = $id=~/^(\S+)/;
           foreach my $p1str (@{$READS1{$id}}) {
             my ($p1,$e1) = split(/\//,$p1str);
             next if $p1 > 0;
             foreach my $p2str (@{$READS2{$id}}) {
               my ($p2,$e2) = split(/\//,$p2str);
               next if $p2 < 0;
               my $size = -$p1-$p2+1;
#               warn "size is $size\n";
               next if ($size > $SIZES{$lib}{'max'} || $size < $SIZES{$lib}{'min'});
               print PAIRS join("\t",1,$id,$p1str,2,$id,$p2str,$size),"\n";
               $MAPPED{$r}++;
#               warn "$r is mapped\n";
             }
          }
      }
   }
   my %LINKREADS;
   foreach my $id (keys %READS1) {
     my ($r,$ctg) = split(/\s/,$id);
     next if defined $MAPPED{$r};
     my ($lib) = $id=~/^(L\d+)/;
     foreach my $p (@{$READS1{$id}}) {
        my ($coord) = $p=~/(\d+)/;
        next if ($coord > $SIZES{$lib}{'max'}); #ignor such cases, since the second in the pair should have also mappped inside given the length of the contig
        push @{$LINKREADS{$r}{'1'}},"$ctg/$p";
        $USED_READS{$r}++;
     }
   }
   foreach my $id (keys %READS2) {
     my ($r,$ctg) = split(/\s/,$id);
     next if defined $MAPPED{$r};
     my ($lib) = $id=~/^(L\d+)/;
     foreach my $p (@{$READS2{$id}}) {
       my ($coord) = $p=~/(\d+)/;
       my $max_possible_size = $CTGSIZE{$ctg} - $coord + 1;
       next if ($max_possible_size > $SIZES{$lib}{'max'}); #ignor such cases, since the second in the pair should have also mappped inside given the length of the contig
       push @{$LINKREADS{$r}{'2'}},"$ctg/$p";
       $USED_READS{$r}++;
     }
   }
   foreach my $r (sort keys %LINKREADS) {
     if (defined $LINKREADS{$r}{'1'}) { print LINKS join("\t","$r/1",@{$LINKREADS{$r}{'1'}}),"\n"; }
     else { print LINKS join("\t","$r/1",""),"\n"; }
     if (defined $LINKREADS{$r}{'2'}) { print LINKS join("\t","$r/2",@{$LINKREADS{$r}{'2'}}),"\n"; }
     else { print LINKS join("\t","$r/2",""),"\n"; }
   }
}

close PAIRS;
close LINKS;

open FQ,">freq";

while (my $f1 = <*.R1.fq>) {
   my $f2=$f1;
   $f2=~s/\.R1\.fq/\.R2\.fq/;
   my $f3=$f1; $f3=~s/\.R1\.fq//;
   open(F1,$f1);
   open(F2,$f2);
   $f1=~s/\.fq$/\.nonmapped\.fq/;
   $f2=~s/\.fq$/\.nonmapped\.fq/;
   open(OUT1,">$f1");
   open(OUT2,">$f2");
   open(OUT3,">$f3.used");
   while (my $h1 = <F1>) {
     my $seq1 = <F1>;
     my $hq1 = <F1>;
     my $qual1 = <F1>;
     my $h2 = <F2>;
     my $seq2 = <F2>;
     my $hq2 = <F2>;
     my $qual2 = <F2>;
     my ($r) = $h1=~/\@(\S+)\//;
     if (defined $FREQ{"$r/1"}) {
       if ($FREQ{"$r/1"} >= $max_freq) {
         print FQ $r,"\t",$FREQ{"$r/1"},"\n";
       }
     }
     elsif (defined $FREQ{"$r/2"}) {
       if ($FREQ{"$r/2"} >= $max_freq) {
         print FQ $r,"\t",$FREQ{"$r/2"},"\n";
       }
     }
     next if defined $MAPPED{$r};
     print OUT1 $h1,$seq1,$hq1,$qual1;
     print OUT2 $h2,$seq2,$hq2,$qual2;
     if (defined $USED_READS{$r}) {
       print OUT3 $h1,$seq1,$hq1,$qual1,$h2,$seq2,$hq2,$qual2;
     }
   }
   close OUT1;
   close OUT2;
   close F1;
   close F2;
}

close FQ;


####################
sub read_multisam {
  my $samfile = shift;
#  warn "parsing $samfile\n";
  open F, $samfile;
  my %READS;
  while (my $bwt = <F>) {
     next unless $bwt=~/AS\:i/;
     my ($read_id,$flag,$ctg,$pos,$mapq,$cigar,$materef,$matepos,$isize,$seq,$qual) = split(/\t/,$bwt);
     $qual=~s/\n//;
     my ($edit_distance) = $bwt=~/NM\:i\:(\d+)/;
     $edit_distance = 100 unless defined $edit_distance;
     next if $edit_distance > $max_distance;
     my $strand = '+';
     if ($flag == 16 || $flag == 272) {
       $strand = '-';
       $pos=$pos+length($seq)-1;
     }
     $FREQ{$read_id}++;
#     warn "freq for $read_id is $FREQ{$read_id}\n";
     $read_id=~s/\/.*//;
     push @{$READS{"$read_id\t$ctg"}}, "$strand$pos/$edit_distance";
  }
  close F;
  return \%READS;
}
##################