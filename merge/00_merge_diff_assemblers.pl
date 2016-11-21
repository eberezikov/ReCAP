#!/usr/bin/perl -w

use strict;


open OUT,">assembled_ctgs.00.fa";
open QUAL,">assembled_ctgs.00.fa.qual";

open(F,"../09_cap_rerun/assembled_ctgs.07.fa") || die "Canpt open cap ctgs\n";
open(Q,"../09_cap_rerun/assembled_ctgs.07.fa.qual") || die "Canpt open cap ctgs\n";
while (my $h = <F>) {
   my $seq = <F>;
   my $hq = <Q>;
   my $qual = <Q>;
   next if length($seq) < 101;
   $h=~s/(>\S+)/$1-cp/;
   $hq=~s/(>\S+)/$1-cp/;
   print OUT $h,$seq;
   print QUAL $hq,$qual;
}
close F;
close Q;
   

open(F,"../09_newbler/assembled_ctgs.07.fa") || die "Canpt open newbler ctgs\n";
open(Q,"../09_newbler/assembled_ctgs.07.fa.qual") || die "Canpt open newbler ctgs\n";
while (my $h = <F>) {
   my $seq = <F>;
   my $hq = <Q>;
   my $qual = <Q>;
   next if length($seq) < 101;
   $h=~s/(>\S+)/$1-nb/;
   $hq=~s/(>\S+)/$1-nb/;
   print OUT $h,$seq;
   print QUAL $hq,$qual;
}
close F;
close Q;





