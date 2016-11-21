#!/usr/bin/perl -w

use strict;

open(FAILED,">failed_ctgs.07.fa");
open(ASS,">assembled_ctgs.07.fa");
open(QUAL,">assembled_ctgs.07.fa.qual");
open(READS,">assembled_ctgs.07.reads");
open(BL2SEQ,">assembled_ctgs.07.bl2seq");

my $n = 0;

while (my $f = <assembly1/*.assembled.fa>) {
   my $data = `cat $f`;
   print ASS $data;
   $data = `cat $f.qual`;
   print QUAL $data;
   $f=~s/\.fa$//;
   $data = `cat $f.reads`;
   print READS $data;
   $data = `cat $f.bl2seq`;
   print BL2SEQ $data;
   $n++;
   warn "$n\n" if $n=~/0$/;
}

while (my $f = <assembly1_re/*.assembled.fa>) {
   my $data = `cat $f`;
   print ASS $data;
   $data = `cat $f.qual`;
   print QUAL $data;
   $f=~s/\.fa$//;
   $data = `cat $f.reads`;
   print READS $data;
   $data = `cat $f.bl2seq`;
   print BL2SEQ $data;
   $n++;
   warn "$n\n" if $n=~/0$/;
}



$n=0;

while (my $f = <assembly1/*.failed.fa>) {
   my $data = `cat $f`;
   print FAILED $data;
}

while (my $f = <assembly1_re/*.failed.fa>) {
   my $data = `cat $f`;
   print FAILED $data;
}



close FAILED;
close ASS;
close QUAL;
close READS;
close BL2SEQ;