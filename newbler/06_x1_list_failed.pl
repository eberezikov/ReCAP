#!/usr/bin/perl -w

use strict;

open OUT,">failed.06";

my $n = 0;
while (my $f = <assembly1/*.o*>) {
   $n++;
   warn "$n: $f\n";
   open(F,$f);
   while(my $l= <F>) {
     if ($l=~/^No assembly/) {
       print OUT "$f\t$l";
     }
   }
}
