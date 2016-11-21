#!/usr/bin/perl -w

use strict;

my $min_reads = 3;

my %COUNTS;
for my $t ('cds','utr5','utr3','nc') {
   open(F,"FINAL/$t.counts") || die "Can't open FINAL/$t.counts\n";
   while (my ($ctg,$r) = split(/\s/,<F>)) {
      $COUNTS{$ctg}{$t}=$r;
   }
}


while (my $f = <FINAL/contigs.cutoff*.list>) {
   warn "$f\n";
   open F,$f;
   open OUT,">$f.classified\n";
   <F>;
   my %TYPES;
   print OUT join("\t",'id','len','type','utr5','cds','utr3','nc'),"\n";
   while (my @l = split(/\s/,<F>)) {
      my $id = $l[0];
      my $len = $l[1];
      $COUNTS{$id}{'cds'}=0 unless defined $COUNTS{$id}{'cds'};
      $COUNTS{$id}{'utr5'}=0 unless defined $COUNTS{$id}{'utr5'};
      $COUNTS{$id}{'utr3'}=0 unless defined $COUNTS{$id}{'utr3'};
      $COUNTS{$id}{'nc'}=0 unless defined $COUNTS{$id}{'nc'};
      my @type;
      if ($COUNTS{$id}{'utr5'} >= $min_reads) {
         push @type, '5UTR';
      }
      if ($COUNTS{$id}{'cds'} >= $min_reads) {
         push @type, 'CDS';
      }
      if ($COUNTS{$id}{'nc'} >= $min_reads && $COUNTS{$id}{'nc'} > 2*$COUNTS{$id}{'cds'}) {
         push @type, 'NC';
      }
      if ($COUNTS{$id}{'utr3'} >= $min_reads) {
         push @type, '3UTR';
      }
      
      push @type, 'LOW' unless @type;

      my $typestr = join('_',@type);
      $TYPES{$typestr}++;

      print OUT join("\t",$id,$len,$typestr,$COUNTS{$id}{'utr5'},$COUNTS{$id}{'cds'},$COUNTS{$id}{'utr3'},$COUNTS{$id}{'nc'}),"\n";
   }
   close F;
   close OUT;
   open OUT,">$f.classified.stats";
   foreach my $type (sort {$TYPES{$b} <=> $TYPES{$a}} keys %TYPES) {
      print OUT $type,"\t",$TYPES{$type},"\n";
   }
   close F;
   close OUT;

}


