#!/usr/bin/perl -w

use strict;

my ($cycledir) = @ARGV;

chdir($cycledir) || die "Can't go to $cycledir\n";

open OUT,">ctgs_out.fa";
open QUAL,">ctgs_out.fa.qual";
open R,">ctgs_out.reads";
open R1,">ctgs_out.assembled";

my $new_ctgs = 0;
my %USED_CTGS;

foreach my $f (<assembly/*.fa>) {
   warn "$f\n";
   open F, $f;
   open Q, "$f.qual";
   while (my $h = <F>) {
     my $seq = <F>;
     my $hq = <Q>;
     my $qual = <Q>;
     print OUT $h,$seq;
     print QUAL $hq,$qual;
     $new_ctgs++;
   }
   close F;
   close Q;
   $f=~s/\.fa$/\.reads/;
   open(F,$f);
   local $/ = "\n\n";
   while (my $readsblock = <F>) {
      $readsblock=~s/(.*)\n//;
      my $readshead = $1;
      my %THISCTGS;
      while ($readsblock=~/(cap\d+\-\S+)/g) {
         $USED_CTGS{$1}++;
         $THISCTGS{$1}++;
      }
      if (scalar keys %THISCTGS > 0) {
        print R $readshead,"\n",$readsblock;
        $readshead=~s/\s.*//;
        print R1 join("\t",$readshead, sort keys %THISCTGS),"\n";
      }
   }
   close F;
   local $/ = "\n";
}

close R;
close R1;


my $old_ctgs = 0;
my $old_nonused_ctgs = 0;
open F,"ctgs_in.fa";
open Q,"ctgs_in.fa.qual";
while (my $h = <F>) {
     my $seq = <F>;
     my $hq = <Q>;
     my $qual = <Q>;
     my ($id) = $h=~/>(\S+)/;
     $old_ctgs++;
     next if defined $USED_CTGS{$id};
     print OUT $h,$seq;
     print QUAL $hq,$qual;
     $old_nonused_ctgs++;
}

close F;
close Q;

close OUT;
close QUAL;

open OUT,">assembled_ctgs.stat";
my $tot_new = $new_ctgs + $old_nonused_ctgs;
my $used = scalar keys %USED_CTGS;
my $tot_old = $used + $old_nonused_ctgs;
print OUT "Total contigs new assembly: $tot_new\n";
print OUT "Total contigs old assembly: $tot_old\n";
print OUT "New contigs: $new_ctgs (",sprintf("%.2f",$new_ctgs/$tot_new*100),"%)\n";
print OUT "Old unchanged contigs: $old_nonused_ctgs (",sprintf("%.2f",$old_nonused_ctgs/$tot_new*100),"% of total, ",sprintf("%.2f",$old_nonused_ctgs/$tot_old*100),"% of old)\n";
print OUT "Old merged contigs: $used (",sprintf("%.2f",$used/$tot_old*100),"% of old)\n";


