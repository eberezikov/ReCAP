#!/usr/bin/perl -w

use strict;

open OUT,">status.05_x0";

while (my $f = <CTG_MERGED1/*.merged>) {
   warn "$f\n";
   open F, $f;
   local $/ = "\n\n";
   while (my $block = <F>) {
      next if $block=~/>/;
      my ($id) = $block=~/^(\S+)/;
      print OUT "$f:\t$id\n";
   }
}
