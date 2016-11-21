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

my $readsdir=$CONF{'RAW_READS'};

chdir 'FINAL';
mkdir 'map';

warn "Running bowtie\n";
`$workdir/run_bowtie3.pl $readsdir $workdir $K $P $NM 2>log.01 1>&2`;

warn "Merging results\n";

`$workdir/merge_mapped2.pl $workdir/FINAL/map 2>log.02 1>&2`;
