#!/usr/bin/perl -w
use strict;
use POSIX qw(tmpnam);
use Config::General qw(ParseConfig);
my %CONF = ParseConfig("CONFIG");

my $RUN_HOST="/mnt/".$CONF{'RUN_HOST'};

my $PRIORITY = $ARGV[0];
   $PRIORITY = 0 unless defined $PRIORITY;
   
my $COVERAGE = $CONF{'DESIRED_COVERAGE'};
   
my $CHUNKS_FOR_MERGING = $CONF{'CHUNKS_FOR_MERGING'};


mkdir "assembly1/failed" unless -e "assembly1/failed";

my %TO_RUN;
open F,"status.06_x1.list";
while (my ($sh,$e) = split(/\s/,<F>)) {
  `mv $e assembly1/failed` if (-e $e);
  $e=~s/\.e/\.o/;
  `mv $e assembly1/failed` if (-e $e);
  $e=~s/.*\///; $e=~s/\..*//;
  $TO_RUN{$sh}=$e;
}


my $workdir = `pwd`;
$workdir=~s/\n//;
$workdir=~s/\/data/$RUN_HOST/;


open  QSUB, ">assembly1/qsub2.sh";
print QSUB "\#!/bin/sh\n\n. /srv/sge/fedot12/common/settings.sh\n\n";

open QSUBDUMMY,">assembly1/qsub2_dummy.sh";
print QSUBDUMMY "\#!/bin/sh\n\n. /srv/sge/fedot12/common/settings.sh\n\n";

my @jobs_to_wait = ();
my $partn=0;
foreach my $p ( keys %TO_RUN ) {
    $partn++;
    my $job_id = $TO_RUN{$p};
    my $pfile = $p;
    warn "submitting $job_id\n";
    print QSUB "qsub -P wga -o $workdir/assembly1 -e $workdir/assembly1 -N $job_id $workdir/assembly1/$pfile.sh\n";
    push @jobs_to_wait, $job_id;
    if (scalar @jobs_to_wait > 100) {
      print QSUBDUMMY "qsub -P wga -o $workdir/assembly1 -e $workdir/assembly1  -sync y -hold_jid ",join(",",@jobs_to_wait)," $workdir/assembly1/dummy_cap.sh\n";
      @jobs_to_wait=();
    }
}
close QSUB;

open SH, ">assembly1/dummy_cap.sh";
print SH "\#!/bin/sh\n\#\$ -S /bin/sh\n";
close SH;

if (@jobs_to_wait) {
 print QSUBDUMMY "qsub -P wga -o $workdir/assembly1 -e $workdir/assembly1  -sync y -hold_jid ",join(",",@jobs_to_wait)," $workdir/assembly1/dummy_cap.sh\n";
}
close QSUBDUMMY;

system("sh $workdir/assembly1/qsub2.sh");
system("sh $workdir/assembly1/qsub2_dummy.sh");


sub get_job_id {
    my $id = tmpnam(); 
    $id =~ s/\/tmp\/file//;
    return $id;
}


