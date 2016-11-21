#!/usr/bin/perl -w

use strict;
use Config::General qw(ParseConfig);
my %CONF = ParseConfig("CONFIG");

my $RUN_HOST="/mnt/".$CONF{'RUN_HOST'};

my $PRIORITY = $ARGV[0];
   $PRIORITY = 0 unless defined $PRIORITY;

my $K = $CONF{'HITS_TO_REPORT'};
my $P = $CONF{'BOWTIE_THREADS'};
my $NM = $CONF{'MAX_MAPPING_DISTANCE'};
my $MAX_LINKED_CTGS = $CONF{'MAX_LINKED_CTGS'};
my $LINKING_THREADS = $CONF{'LINKING_THREADS'};
my $CHUNKS_FOR_MERGING = $CONF{'CHUNKS_FOR_MERGING'};
my $ASSEMBLY2_CHUNKS = $CONF{'ASSEMBLY2_CHUNKS'};


my $workdir = `pwd`;
$workdir=~s/\n//;
$workdir=~s/\/data/$RUN_HOST/;

my $cycle;

if (-e 'MERGING') {
  if (-e 'MERGING/CYCLE_TO_RUN') {
     $cycle = `cat MERGING/CYCLE_TO_RUN`;
     $cycle=~s/\n//;
  }
  else {
    die "No cycle info found\n\n";
  }
}
else {
  die "No MERGING folde exists\n";
}

if ($cycle=~/FINAL/) {
  print "Merging has been already complete.\n";
  exit;
}

if ($cycle == 1) {
  print "This script can be run only from the second cycle. Current cycle us $cycle\n";
  exit;
}

my %LIBNAME;
open(F,'LIBNAMES') || die "Can't read LIBNAMES\n";;
while (my @l = split(/\s/,<F>)) {
   $LIBNAME{$l[1]}=$l[0];
   $LIBNAME{$l[0]}=$l[1];
}
close F;

chdir 'MERGING';


my $prev = $cycle - 1;

my %FASTQ;
for my $c (1..$prev) {
 while (my $f = <c$c/map/*/*.used_in_links>) {
    read_fastq_used_in_links($f);
 }
}
print "Done reading initial fastq\n";


my %CTG;
for my $c (1..$prev) {
  read_contigs("c$c/ctgs_in.fa");
}

my %MAPPING;
for my $c (1..$prev) {
 while (my $lnk = <c$c/map/*/*.mapped_links>) {
    read_mapping($lnk);
 }
}
print "Done reading initial mapping\n";

my $trimming = 0;
my $trimming_switch = 0;

while (1) {
   last if ($cycle eq 'FINAL' && $trimming == 1);
   $prev = $cycle - 1;
   mkdir "c$cycle";
   mkdir "c$cycle/map";
   if ($trimming_switch == 1) {
     `cp c$prev/ctgs_out.nr.fa c$cycle/ctgs_in.fa`;
     `cp c$prev/ctgs_out.nr.fa.qual c$cycle/ctgs_in.fa.qual`;
     `cp c$prev/pairs.nofrequent c$cycle/pairs.nofrequent`;
     chdir "c$cycle";
     `ln -s pairs.nofrequent pairs.nofrequent.new`;
     $trimming_switch = 0;
   }
   else {
     `cp c$prev/ctgs_out.nr.new.fa c$cycle/ctgs_in.fa`;
     `cp c$prev/ctgs_out.nr.new.fa.qual c$cycle/ctgs_in.fa.qual`;
     chdir "c$cycle";
     print "Building bowtie index for cycle $cycle\n";
     `/usr/local/bowtie/2.2.4/bowtie2-build ctgs_in.fa ctgs_in.fa 2>log.00 1>&2`;
   
     if ($cycle == 2) {
        print "Splitting reads\n";
       `$workdir/split_reads.pl $workdir/MERGING/c1/map $workdir/MERGING/c2/map 2>log.01 1>&2`;
     }
   
     print "Running read mapping\n";
     `$workdir/run_bowtie2.pl $workdir/MERGING/c2/map $workdir c$cycle $K $P $NM 2>log.02 1>&2`;

     print "Merging results\n";
     `$workdir/merge_mapped.pl $workdir/MERGING/c$cycle/map 2>log.03 1>&2`;

     print "Calculating linked ctg pairs\n";
     `$workdir/calculate_linked_pairs_recursive.pl $workdir/MERGING/c$cycle $workdir/LIBNAMES $MAX_LINKED_CTGS 2>log.04 1>&2`;

     print "Filtering frequent contigs\n";
     `$workdir/filter_frequent_contigs.pl $workdir/MERGING/c$cycle/pairs $workdir/MERGING/c$cycle/pairs.nofrequent 5  2>log.05 1>&2`;
     `$workdir/get_new_pairs_recursive.pl $workdir/MERGING/c$cycle 2>log.06 1>&2`;
   }


    while (my $f = <map/*/*.used_in_links>) {
      read_fastq_used_in_links($f);
    }

    read_contigs('ctgs_in.fa');

    while (my $f = <map/*/*.mapped_links>) {
      read_mapping($f);
    }

    make_ctg_reads();

    print "Running read merging\n";
    `$workdir/run_read_merging.pl $workdir $workdir/MERGING/c$cycle $CHUNKS_FOR_MERGING 2>log.08 1>&2`;

    print "Running ctg assembly\n";
    `$workdir/run_cap_est_assemblies2_pairs.pl $workdir $cycle $trimming $ASSEMBLY2_CHUNKS 2>log.09 1>&2`;

    print "Merging ctg assembly\n";
    `$workdir/merge_est_results.pl $workdir/MERGING/c$cycle 2>log.10 1>&2`;

    print "Merging assemblies\n";
    `$workdir/merge_assemblies.pl $cycle 2>log.11 1>&2`;

    print "Removing redundancy\n";
    `/home/fedot/src/cd-hit-v4.6.1-2012-08-27/cd-hit-est -i ctgs_out.all.fa -o ctgs_out.nr.fa -r 0 -c 1 -T $LINKING_THREADS -M 0 2>log.12 1>&2`;

    print "Adding quality data\n";
    `$workdir/add_quality_to_nr.pl $workdir $cycle 2>log.13 1>&2`;

    print "Calculating statistics\n";
    `$workdir/make_n50.pl $workdir/MERGING/c$cycle/ctgs_out.nr.fa >$workdir/MERGING/c$cycle/ctgs_out.nr.fa.n50 2>log.14`;

    my $newcycle = $cycle + 1;
    print "Extracting new contigs\n";
    `$workdir/get_new_contigs.pl $newcycle  $workdir/MERGING/c$cycle/ctgs_out.nr.fa $workdir/MERGING/c$cycle/ctgs_out.nr.new.fa 2>log.15`;

    print "All done for cycle $cycle\n";
    
    chdir "..";
    $cycle = `cat CYCLE_TO_RUN`;
    $cycle=~s/\n//g;
    if ($cycle=~/FINAL/) {
      if ($trimming == 0) {
        `cp CYCLE_TO_RUN CYCLES_NOT_TRIMMED`;
        $trimming = 1;
        $trimming_switch = 1;
        $cycle=~s/FINAL//;
      }
      else {
        $cycle='FINAL' 
      }
    }
#    last;
}


###############

sub make_ctg_reads {
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
  print "There are $totctg contigs to consider\n";
  print "But ", scalar keys %CTG_TO_SKIP, " were already generated, skipping them for the merging step\n";

  foreach my $ctg (keys %CTG_TO_SKIP) {
    delete $CTG_TO_CONSIDER{$ctg};
  }

  mkdir 'CTG_SEQ';

  print "Writing reads for contigs\n";

  my $current = 0;
  foreach my $c (keys %CTG_TO_CONSIDER) {
     $current++;
     print "$current/$totctg.\n";# if $current=~/00$/;
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
}

##################

sub read_fastq_used_in_links {
    my $f = shift;
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

###########

sub read_contigs {
  my $f = shift;
  warn "Reading contigs from $f\n";
  open(F,$f) || die "Can't open $f\n";
  open(Q,"$f.qual") || die "Can't open $f.qual\n";
  while (my $h = <F>) {
    my $seq = <F>;
    my $hq = <Q>;
    my $qual = <Q>;
    my ($id) = $h=~/>(\S+)/;
    $CTG{$id}=$seq.$qual;;
  }
  close F;
  close Q;
}

##########

sub read_mapping {
   my $lnk = shift;
   warn "Reading mapping info from $lnk\n";
   open(F,$lnk);
   my $totmap = `wc -l $lnk`;
      $totmap=~s/\s.*\n//;
   my $curr_map = 0;
   while (my ($r,@ctgs) = split(/\s/,<F>)) {
        $r=~s/\/.*//;
        foreach my $c (@ctgs) {
          $c=~s/\/.*//;
          push @{$MAPPING{$c}},$r;
        }
        $curr_map++;
        warn "$lnk: $curr_map/$totmap ",sprintf("%.2f",$curr_map/$totmap*100),"%\r\n" if ($curr_map=~/00000$/);
   }
   close F;
}

#########
