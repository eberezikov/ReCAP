#!/usr/bin/perl -w
use strict;
use POSIX qw(tmpnam);
use Config::General qw(ParseConfig);
my %CONF = ParseConfig("CONFIG");

my $RUN_HOST="/mnt/".$CONF{'RUN_HOST'};

my $PRIORITY = $ARGV[0];
   $PRIORITY = 0 unless defined $PRIORITY;
   
my $MERGING_MAX_READS = $CONF{'MERGING_MAX_READS_RE'};
my $MERGING_MIN_READS = $CONF{'MERGING_MIN_READS_RE'};
my $MERGING_WINDOW = $CONF{'MERGING_WINDOW_RE'};
   
my $CHUNKS_FOR_MERGING = $CONF{'CHUNKS_FOR_MERGING'};

if (-e 'CTG_MERGED1_RE') {
   die "CTG_MERGED1_RE already exists!\n";
}
else {
  mkdir 'CTG_MERGED1_RE';
}

my @to_run;
while (my $f = <ctgs_rerun/*.ctg>) {
  push @to_run, $f;
}

my $workdir = `pwd`;
$workdir=~s/\n//;
$workdir=~s/\/data/$RUN_HOST/;


open  QSUB, ">CTG_MERGED1_RE/qsub1.sh";
print QSUB "\#!/bin/sh\n\n. /srv/sge/fedot12/common/settings.sh\n\n";

open QSUBDUMMY,">CTG_MERGED1_RE/qsub1_dummy.sh";
print QSUBDUMMY "\#!/bin/sh\n\n. /srv/sge/fedot12/common/settings.sh\n\n";

my @jobs_to_wait = ();
my $partn=0;
foreach my $p ( @to_run ) {
    $partn++;
    my $job_id = join( '_', "mrg",$partn,get_job_id() );
    my $pfile = "mrg\_$partn";
    warn "submitting $job_id\n";
    open  SH, ">CTG_MERGED1_RE/$pfile.sh";
    print SH "\#!/bin/sh\n\#\$ -S /bin/sh\n";
    print SH "mkdir /tmp/$job_id\n";
    print SH "cd /tmp/$job_id\n";
    print SH "uname -n >$pfile.host\n";
    print SH "date >>$pfile.host\n";
    print SH "cp $workdir/LIBNAMES ./\n";
    print SH "cp $workdir/$p ./\n";
    $p=~s/.*\///;
    print SH "$workdir/normalize_coverage_and_merge_reads_multi.pl $p LIBNAMES $MERGING_MAX_READS $MERGING_WINDOW\n";
    print SH "mv *.merged *.reduction $workdir/CTG_MERGED1_RE\n";
    print SH "date >>$pfile.host\n";
    print SH "mv $pfile.host $workdir/CTG_MERGED1_RE\n";
    print SH "cd ..\n";
    print SH "rm -R $job_id\n"; 
    close SH;
    print QSUB "qsub -P wga -o $workdir/CTG_MERGED1_RE -e $workdir/CTG_MERGED1_RE -N $job_id $workdir/CTG_MERGED1_RE/$pfile.sh\n";
    push @jobs_to_wait, $job_id;
    if (scalar @jobs_to_wait > 100) {
      print QSUBDUMMY "qsub -P wga -o $workdir/CTG_MERGED1_RE -e $workdir/CTG_MERGED1_RE  -sync y -hold_jid ",join(",",@jobs_to_wait)," $workdir/CTG_MERGED1_RE/dummy_cap.sh\n";
      @jobs_to_wait=();
    }
}
close QSUB;

open SH, ">CTG_MERGED1_RE/dummy_cap.sh";
print SH "\#!/bin/sh\n\#\$ -S /bin/sh\n";
close SH;

if (@jobs_to_wait) {
 print QSUBDUMMY "qsub -P wga -o $workdir/CTG_MERGED1_RE -e $workdir/CTG_MERGED1_RE  -sync y -hold_jid ",join(",",@jobs_to_wait)," $workdir/CTG_MERGED1_RE/dummy_cap.sh\n";
}
close QSUBDUMMY;

system("sh $workdir/CTG_MERGED1_RE/qsub1.sh");
system("sh $workdir/CTG_MERGED1_RE/qsub1_dummy.sh");


sub get_job_id {
    my $id = tmpnam(); 
    $id =~ s/\/tmp\/file//;
    return $id;
}

