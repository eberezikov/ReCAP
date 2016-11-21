#!/usr/bin/perl -w

use strict;

my ($READSDIR) = @ARGV;


while (my $dir  = <$READSDIR/*>) {
  my ($lib) = $dir=~/$READSDIR\/(\S+)/;

  open OUT1,">$dir/$lib.nonmapped.R1.fastq";
  open OUT2,">$dir/$lib.nonmapped.R2.fastq";
  while (my $f1 = <$dir/*.R1.nonmapped.fq>) {
       warn "$f1\n";
       my $seq1 = `cat $f1`;
       print OUT1 $seq1;
       my $f2=$f1;
       $f2=~s/\.R1\./\.R2\./;
       my $seq2 = `cat $f2`;
       print OUT2 $seq2;
  }
  close OUT1;
  close OUT2;

  open OUT,">$dir/$lib.mapped_reads";
  while (my $f = <$dir/*.map>) {
    my $data = `cat $f`;
    print OUT $data;
    `rm $f`;
  }
  close OUT;

  open OUT,">$dir/$lib.mapped_links";
  while (my $f = <$dir/*.lnk>) {
    my $data = `cat $f`;
    print OUT $data;
    `rm $f`;
  }
  close OUT;

  open OUT,">$dir/$lib.frequent";
  while (my $f = <$dir/*.frq>) {
    my $data = `cat $f`;
    print OUT $data;
    `rm $f`;
  }
  close OUT;

  open OUT,">$dir/$lib.used_in_links";
  while (my $f = <$dir/*.used>) {
    my $data = `cat $f`;
    print OUT $data;
    `rm $f`;
  }
  close OUT;

  while (my $f = <$dir/*.e*>) { `rm $f`; }
  while (my $f = <$dir/*.o*>) { `rm $f`; }
  while (my $f = <$dir/*.pe*>) { `rm $f`; }
  while (my $f = <$dir/*.po*>) { `rm $f`; }
  while (my $f = <$dir/*.sh>) { `rm $f`; }
  while (my $f = <$dir/*.host>) { `rm $f`; }
  while (my $f = <$dir/*.nonmapped.fq>) { `rm $f`; }

}
