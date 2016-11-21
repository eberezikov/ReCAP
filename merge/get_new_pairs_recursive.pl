#!/usr/bin/perl -w

use strict;

my ($cycledir) = @ARGV;

$cycledir=~s/(\d+)$//;
my $last_cycle = $1;


my %OLD;

for my $c (1..$last_cycle) {
  if (-e "$cycledir$c/pairs.nofrequent.new") {
     warn "readng previous pairs from $cycledir$c/pairs.nofrequent.new\n";
     open(F,"$cycledir$c/pairs.nofrequent.new");
     while (my @g = split(/\s/,<F>)) {
       $OLD{"$g[0] $g[1]"}++;
     }
     close F;
  }
}

warn "writing down $cycledir$last_cycle/pairs.nofrequent.new\n";
open(F,"$cycledir$last_cycle/pairs.nofrequent");
open(OUT,">$cycledir$last_cycle/pairs.nofrequent.new");
while (my @g = split(/\s/,<F>)) {
  next if defined $OLD{"$g[0] $g[1]"};
  print OUT join("\t",@g),"\n";
}
close F;
close OUT;
