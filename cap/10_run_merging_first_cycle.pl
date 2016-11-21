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
my $readsdir;


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
  warn "Initializing first merging cycle\n";
  $cycle = 1;
  mkdir 'MERGING';
  mkdir 'MERGING/CTG_MERGED';
}

if ($cycle != 1) {
  die "This script is only for the first cycle!\n";
}

if ($cycle=~/FINAL/) {
  warn "Merging complete. Generating final output\n";

}
else {
  warn "Generating script to run cycle $cycle\n";
  mkdir "MERGING/c$cycle";
  mkdir "MERGING/c$cycle/map";
  open(SH,">MERGING/c$cycle/run.sh");
  print SH "cd $workdir/MERGING/c$cycle\n";
  print SH "echo \"Building bowtie index for cycle $cycle\"\n";
  print SH "/usr/local/bowtie/2.2.4/bowtie2-build ctgs_in.fa ctgs_in.fa 2>log.00 1>&2\n\n";
  if ($cycle == 1) {
    mkdir 'MERGING/CTG_MERGED' unless -e 'MERGING/CTG_MERGED';
   `cp assembled_ctgs.07.nr99.fa MERGING/c$cycle/ctgs_in.fa`;
   `cp assembled_ctgs.07.nr99.fa.qual MERGING/c$cycle/ctgs_in.fa.qual`;
   $readsdir=$CONF{'RAW_READS'};
   print SH "echo \"Running read mapping\"\n";
   print SH "$workdir/run_bowtie2.pl $readsdir $workdir c$cycle $K $P $NM 2>log.02 1>&2\n\n";
  }
  else {
    my $prev = $cycle - 1;
    `cp MERGING/c$prev/ctgs_out.nr.new.fa MERGING/c$cycle/ctgs_in.fa`;
    `cp MERGING/c$prev/ctgs_out.nr.new.fa.qual MERGING/c$cycle/ctgs_in.fa.qual`;
    $readsdir = "$workdir/MERGING/c$prev/map";
    print SH "echo \"Splitting reads\"\n";
    print SH "$workdir/split_reads.pl $readsdir $workdir/MERGING/c$cycle/map 2>log.01 1>&2\n\n";
    print SH "echo \"Running read mapping\"\n";
    print SH "$workdir/run_bowtie2.pl $workdir/MERGING/c$cycle/map $workdir c$cycle $K $P $NM 2>log.02 1>&2\n\n";
  }

  print SH "echo \"Merging results\"\n";
  print SH "$workdir/merge_mapped.pl $workdir/MERGING/c$cycle/map 2>log.03 1>&2\n\n";

  print SH "echo \"Calculating linked ctg pairs\"\n";
  print SH "$workdir/calculate_linked_pairs_recursive.pl $workdir/MERGING/c$cycle $workdir/LIBNAMES $MAX_LINKED_CTGS 2>log.04 1>&2\n\n";

  print SH "echo \"Filtering frequent contigs\"\n";
  print SH "$workdir/filter_frequent_contigs.pl $workdir/MERGING/c$cycle/pairs $workdir/MERGING/c$cycle/pairs.nofrequent 5  2>log.05 1>&2\n\n";


  if ($cycle == 1) {
    print SH "ln -s $workdir/MERGING/c$cycle/pairs.nofrequent $workdir/MERGING/c$cycle/pairs.nofrequent.new\n\n";
  }
  else {
    my $prev = $cycle - 1;
    print SH "$workdir/get_new_pairs_recursive.pl $workdir/MERGING/c$cycle 2>log.06 1>&2\n\n";
  }

  print SH "echo \"Generating ctg seqs\"\n";
  print SH "$workdir/make_ctg_reads_in_memory_recursive.pl $workdir/MERGING/c$cycle $workdir/LIBNAMES 2>log.07 1>&2\n\n";


  print SH "echo \"Running read merging\"\n";
  print SH "$workdir/run_read_merging.pl $workdir $workdir/MERGING/c$cycle $CHUNKS_FOR_MERGING 2>log.08 1>&2\n\n";

  print SH "echo \"Running ctg assembly\"\n";
  print SH "$workdir/run_cap_est_assemblies2_pairs.pl $workdir $cycle $ASSEMBLY2_CHUNKS 2>log.09 1>&2\n\n";

  print SH "echo \"Merging ctg assembly\"\n";
  print SH "$workdir/merge_est_results.pl $workdir/MERGING/c$cycle 2>log.10 1>&2\n\n";

  print SH "echo \"Removing redundancy \"\n";
  print SH "/home/fedot/src/cd-hit-v4.6.1-2012-08-27/cd-hit-est -i ctgs_out.fa -o ctgs_out.nr.fa -r 0 -c 1 -T $LINKING_THREADS -M 0 2>log.11 1>&2\n\n";

  if ($cycle == 1) {
    print SH "ln -s $workdir/MERGING/c$cycle/ctgs_out.fa $workdir/MERGING/c$cycle/ctgs_out.all.fa\n\n";
    print SH "ln -s $workdir/MERGING/c$cycle/ctgs_out.fa.qual $workdir/MERGING/c$cycle/ctgs_out.all.fa.qual\n\n";
  }


  print SH "echo \"Adding quality data \"\n";
  print SH "$workdir/add_quality_to_nr.pl $workdir $cycle 2>log.12 1>&2\n\n";

  print SH "echo \"Calculating statistics \"\n";
  print SH "$workdir/make_n50.pl $workdir/MERGING/c$cycle/ctgs_out.nr.fa >$workdir/MERGING/c$cycle/ctgs_out.nr.fa.n50 2>log.13\n\n";

  my $newcycle = $cycle + 1;
  print SH "echo \"Extracting new contigs \"\n";
  print SH "$workdir/get_new_contigs.pl $newcycle  $workdir/MERGING/c$cycle/ctgs_out.nr.fa $workdir/MERGING/c$cycle/ctgs_out.nr.new.fa 2>log.14\n\n";


  print SH "echo \"All done\"\n";
  close SH;
  
  warn "\n\nNow run MERGING/c$cycle/run.sh\n\n";


}