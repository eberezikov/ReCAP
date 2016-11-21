#!/usr/bin/perl -w

use strict;

my ($cycledir,$libnames) = @ARGV;

$cycledir=~s/(\d+)$//;
my $last_cycle = $1;



my %LIBNAME;
open(F,$libnames) || die "Can't read $libnames\n";;
while (my @l = split(/\s/,<F>)) {
   $LIBNAME{$l[1]}=$l[0];
   $LIBNAME{$l[0]}=$l[1];
}
close F;

chdir("$cycledir$last_cycle") || die "Can't go to $cycledir$last_cycle\n";

my %CTG_TO_CONSIDER;
open(F,"pairs.nofrequent.new") || die "Can't open groups\n";
while (my (@ctgs) = split(/\s/,<F>)) {
    $CTG_TO_CONSIDER{$ctgs[0]}++;
    $CTG_TO_CONSIDER{$ctgs[1]}++;
}
close F;

my %CTG_TO_SKIP;
#skipping previously generates ctgs\n";
foreach my $ctg (keys %CTG_TO_CONSIDER) {
  if (-e "../CTG_MERGED/$ctg.ctg.merged") {
    $CTG_TO_SKIP{$ctg}++;
  }
}

my $totctg = scalar keys %CTG_TO_CONSIDER;
warn "There are $totctg contigs to consider\n";
warn "But ", scalar keys %CTG_TO_SKIP, " were already generated, skipping them for the merging step\n";

foreach my $ctg (keys %CTG_TO_SKIP) {
  delete $CTG_TO_CONSIDER{$ctg};
}

my %CTG;
for my $c (1..$last_cycle) {
  warn "Reading contigs from $cycledir$c/ctgs_in.fa\n";
  open(F,"$cycledir$c/ctgs_in.fa") || die "Can't open $cycledir$c/ctgs_in.fa\n";
  open(Q,"$cycledir$c/ctgs_in.fa.qual") || die "Can't open $cycledir$c/ctgs_in.fa.qual\n";
  while (my $h = <F>) {
    my $seq = <F>;
    my $hq = <Q>;
    my $qual = <Q>;
    my ($id,$len) = $h=~/>(\S+) (\d+)/;
    $CTG{$id}=$seq.$qual if (defined $CTG_TO_CONSIDER{$id});
  }
  close F;
  close Q;
}

my %MAPPING;
for my $c (1..$last_cycle) {
  while (my $lnk = <$cycledir$c/map/*/*.mapped_links>) {
    warn "Reading mapping info from $lnk\n";
    open(F,$lnk);
    my $totmap = `wc -l $lnk`;
    $totmap=~s/\s.*\n//;
    my $curr_map = 0;
    while (my ($r,@ctgs) = split(/\s/,<F>)) {
      $r=~s/\/.*//;
      foreach my $c (@ctgs) {
        $c=~s/\/.*//;
        push @{$MAPPING{$c}},$r if (defined $CTG_TO_CONSIDER{$c});
      }
      $curr_map++;
      warn "$lnk: $curr_map/$totmap ",sprintf("%.2f",$curr_map/$totmap*100),"%\r\n" if ($curr_map=~/00000$/);
    }
  }
}
warn "done reading mapping info\n";
close F;

my %FASTQ;
for my $c (1..$last_cycle) {
  while (my $f = <$cycledir$c/map/*/*.used_in_links>) {
    warn "Reading fastq from $f\n";
    my ($libname) = $f=~/map\/(\S+)\//;
    my $libid = $LIBNAME{$libname};
    my $totreads =  `wc -l $f`;
    $totreads=~s/\s.*\n//;
    $totreads=$totreads/8;
    my $curr_read = 0;
    open(F,$f);
    while (my $h1 = <F>) {
       my $seq1=<F>;
       my $hq1=<F>;
       my $qual1=<F>;
       my $h2 = <F>;
       my $seq2=<F>;
       my $hq2=<F>;
       my $qual2=<F>;
       my ($id) = $h1=~/r(\d+)/;
       $FASTQ{$libid}->[$id]=$h1.$seq1.$hq1.$qual1.$h2.$seq2.$hq2.$qual2;
       $curr_read++;
       warn "$f: $curr_read/$totreads ",sprintf("%.2f",$curr_read/$totreads*100),"%\n" if ($curr_read=~/00000$/);
    }
    close F;
  }
}

warn "Done reading fastq\n";

mkdir 'CTG_SEQ';

warn "Writing reads for contigs\n";

my $current = 0;
foreach my $c (keys %CTG_TO_CONSIDER) {
   $current++;
   warn "$current/$totctg.\n";# if $current=~/00$/;
   open OUT,">CTG_SEQ/$c.ctg";
   my %READS;
   print OUT $c,"\n",$CTG{$c};
   if (defined $MAPPING{$c}) {
      foreach my $r (@{$MAPPING{$c}}) {
        $READS{$r}++;
      }
   }
   foreach my $r (keys %READS) {
      my ($lib,$nreads) = $r=~/(L\d+)r(\d+)/;
      print OUT $FASTQ{$lib}->[$nreads];
   }
   close OUT;
}

