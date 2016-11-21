#!/usr/bin/perl -w

use strict;

my ($cycledir,$libnames,$MAX_LINKED_CTGS) = @ARGV;

$MAX_LINKED_CTGS = 1 unless defined $MAX_LINKED_CTGS;

my %SIZES;
open(F,$libnames) || die "Can't open libnames\n";;
while (my @l = split(/\s/,<F>)) {
   $SIZES{$l[1]}{'min'}=$l[2];
   $SIZES{$l[1]}{'max'}=$l[3];
   $SIZES{$l[1]}{'origname'}=$l[0];
}
close F;

$cycledir=~s/(\d+)$//;
my $last_cycle = $1;


my %CTGLEN;
for my $c (1..$last_cycle) {
  warn "reading contig lengths from $cycledir$c/ctgs_in.fa\n";
  open(F,"$cycledir$c/ctgs_in.fa") || die "Can't open input contigs\n";
  while (my $h = <F>) {
    <F>;
    my ($id,$len) = $h=~/>(\S+) (\d+)/;
    $CTGLEN{$id}=$len;
  }
  close F;
}

my %ASSEMBLED_CTGS;
for my $c (1..$last_cycle) {
  if (-e "$cycledir$c/ctgs_out.assembled") {
    warn "reading assembled ctgs from $cycledir$c/ctgs_out.assembled\n";
    open(F,"$cycledir$c/ctgs_out.assembled");
    while (my ($new,@old) = split(/\s/,<F>)) {
       foreach my $id (@old) {
         $ASSEMBLED_CTGS{$id}++;
       }
    }
  }
}


my %PAIRS;
for my $c (1..$last_cycle) {
 foreach my $f (<$cycledir$c/map/*/*.mapped_links>) {
   warn "parsing links from $f\n";
   open(F,$f);
   while (my ($id1,@links1) = split(/\s/,<F>)) {
     my ($id2,@links2) = split(/\s/,<F>);
     next unless @links1;
     next unless @links2;
     next if (scalar @links1 > $MAX_LINKED_CTGS || scalar @links2 > $MAX_LINKED_CTGS); #drop contigs that connect to many other contigs
     foreach my $l1 (@links1) {
       my ($c1,$pos1,$m1) = split(/\//,$l1);
       next if defined $ASSEMBLED_CTGS{$c1};
       next if $pos1 > 0; # minus strand only;
       next if $m1 > 3; # max mismatches
       my ($lib) = $id1=~/(L\d+)/;
       foreach my $l2 (@links2) {
          my ($c2,$pos2,$m2) = split(/\//,$l2);
          next if defined $ASSEMBLED_CTGS{$c2};
          next if $pos2 < 0; # plus strand only;
          next if $m2 > 3; # max mismatches
          my $len2 = $CTGLEN{$c2} - $pos2 + 1;
          my $len1 = -$pos1;
          my ($max_overlap) = $len2 + $len1;
          next if ($max_overlap  < $SIZES{$lib}{'min'});
          my ($id) = $id1=~/(\S+)\//;
          push @{$PAIRS{"$c2 $c1"}},$id;
       }
     }
   }
   close F;
 }
}

#my %USED_CTGS;

warn "writing overlapping pairs\n";
open(OUT1,">$cycledir$last_cycle/pairs") || die "Can't open pairs for writing\n";;
open(OUT2,">$cycledir$last_cycle/pairs.same") || die "Can't open pairs.same for writing\n";
foreach my $p (sort {scalar @{$PAIRS{$b}} <=> scalar @{$PAIRS{$a}}} keys %PAIRS) {
   my ($p1,$p2) = split(/\s/,$p);
   if ($p1 ne $p2) {
#     next if defined $USED_CTGS{$p1};
#     next if defined $USED_CTGS{$p2};
     print OUT1 join(" ",$p,scalar @{$PAIRS{$p}}),"\n";
#     $USED_CTGS{$p1}++;
#     $USED_CTGS{$p2}++;
   }
   else {
#     print OUT2 join(" ",$p,scalar @{$PAIRS{$p}},@{$PAIRS{$p}}),"\n";
     print OUT2 join(" ",$p,scalar @{$PAIRS{$p}}),"\n";
   }
}
