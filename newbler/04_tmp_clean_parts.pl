#!/usr/bin/perl -w

use strict;

for my $lib (<LIBS_LS/*>) {
  warn "$lib\n";
  `rm $lib/*.ctg`;
}


#for my $p (5000..10000) {
#  warn "$p\n";
#  `rm LIBS_LS/*/*.part$p.ctg`;
#}
