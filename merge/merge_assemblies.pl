#!/usr/bin/perl -w

use strict;

my ($last_cycle) = @ARGV;

my %TO_SKIP;

for my $c (1..$last_cycle) {
  open(F,"../c$c/ctgs_out.assembled") || die "Can't open ../c$c/ctgs_out.assembled\n";
  warn "reading assembled ctgs from ../c$c/ctgs_out.assembled \n";
  while (my ($new,@old) = split(/\s/,<F>)) {
     foreach my $id (@old) {
       $TO_SKIP{$id}++;
     }
  }
  close F;
}

my %PRINTED;

open FA,">ctgs_out.all.fa";
open QV,">ctgs_out.all.fa.qual";
for my $c (1..$last_cycle-1) {
  warn "parsing ../c$c/ctgs_out.nr.fa\n";
  open(F,"../c$c/ctgs_out.nr.fa") || die "Can't open ../c$c/ctgs_out.nr.fa\n";
  open(Q,"../c$c/ctgs_out.nr.fa.qual") || die "Can't open ../c$c/ctgs_out.nr.fa.qual\n";
  while (my $h = <F>) {
    my $seq = <F>;
    my $hq = <Q>;
    my $qual = <Q>;
    my ($id) = $h=~/>(\S+)/;
    next if defined $TO_SKIP{$id};
    next if defined $PRINTED{$id};
    print FA $h,$seq;
    print QV $hq,$qual;
    $PRINTED{$id}++;
  }
  close F;
  close Q;
}

open(F,"ctgs_out.fa") || die "Can't open ctgs_out.fa\n";
open(Q,"ctgs_out.fa.qual") || die "Can't open ctgs_out.fa.qual\n";
warn "parsing latest ctgs_out.fa\n";
while (my $h = <F>) {
    my $seq = <F>;
    my $hq = <Q>;
    my $qual= <Q>;
    my ($id) = $h=~/>(\S+)/;
    next if defined $TO_SKIP{$id};
    next if defined $PRINTED{$id};
    print FA $h,$seq;
    print QV $hq,$qual;
    $PRINTED{$id}++;
}
close F;
close Q;

close FA;
close QV;

