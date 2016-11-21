#!/usr/bin/perl -w
use strict;
use POSIX qw(tmpnam);
use Config::General qw(ParseConfig);
my %CONF = ParseConfig("CONFIG");

my $RUN_HOST="/mnt/".$CONF{'RUN_HOST'};

my $PRIORITY = $ARGV[0];
   $PRIORITY = 0 unless defined $PRIORITY;

my $GENOME = $CONF{'BWA_DB'};

my $K = $CONF{'HITS_TO_REPORT'};
my $P = $CONF{'BOWTIE_THREADS'};
my $NM = $CONF{'MAX_EDIT'};
my $ctgparts = $CONF{'CTGPARTS'};


my $workdir = `pwd`;
$workdir=~s/\n//;
$workdir=~s/\/data/$RUN_HOST/;

my $readsdir=$workdir;
$readsdir=~s/(.*\/).*/${1}00_RAW/;


mkdir "LIBS_LS";
while (my $d = <$readsdir/LIBS_LS/*>) {
  $d=~s/.*\///;
  mkdir "LIBS_LS/$d";
}

my @to_run;
while (my $f = <$readsdir/LIBS_LS/*/*.R1.fq>) {
  push @to_run, $f;
#  last;
}

open  QSUB, ">qsub2.sh";
print QSUB "\#!/bin/sh\n\n. /srv/sge/fedot12/common/settings.sh\n\n";


my @jobs_to_wait = ();
my $partn=0;
foreach my $p ( @to_run ) {
    $partn++;
    my $job_id = join( '_', "bwt",$partn,get_job_id() );
    my $pfile = $p; $pfile=~s/.*\///; $pfile=~s/\.R1\.fq//;
    my ($pdir) = $p=~/(LIBS_LS\/\S+)\//;
    warn "submitting $pfile\n";
    open  SH, ">$pdir/$job_id.sh";
    print SH "\#!/bin/sh\n\#\$ -S /bin/sh\n";
    print SH "mkdir /tmp/$job_id\n";
    print SH "cd /tmp/$job_id\n";
    print SH "uname -n >$pfile.host\n";
    print SH "cp $readsdir/$pdir/$pfile.R1.fq $readsdir/$pdir/$pfile.R2.fq ./\n";
    print SH "date >>$pfile.host\n";
    print SH "$workdir/run_bowtie_global_filter_best_and_distance_and_strand.pl $GENOME $pfile.R1.fq $pfile.R1.sam $pfile.R1.err $P $K $NM R\n";
    print SH "$workdir/run_bowtie_global_filter_best_and_distance_and_strand.pl $GENOME $pfile.R2.fq $pfile.R2.sam $pfile.R2.err $P $K $NM F\n";
    print SH "$workdir/sam_to_parts.pl $pfile.R1.sam $pfile.R2.sam $pfile.R1.fq $pfile.R2.fq $pfile.sams $ctgparts $K\n";
    print SH "date >>$pfile.host\n";
    print SH "cp $pfile.sams $pfile.host $pfile.*.err $workdir/$pdir\n";
    print SH "cd ..\n";
    print SH "rm -R /tmp/$job_id\n"; 
    close SH;
    print QSUB "qsub -P wga -pe threaded $P -o $workdir/$pdir -e $workdir/$pdir -N $job_id $workdir/$pdir/$job_id.sh\n";
    push @jobs_to_wait, $job_id;
}

open SH, ">dummy_bwt.sh";
print SH "\#!/bin/sh\n\#\$ -S /bin/sh\n";
close SH;

print QSUB "qsub -P wga -o $workdir -e $workdir  -sync y -hold_jid ",join(",",@jobs_to_wait)," $workdir/dummy_bwt.sh";
close QSUB;

system("sh $workdir/qsub2.sh");

sub get_job_id {
    my $id = tmpnam(); 
    $id =~ s/\/tmp\/file//;
    return $id;
}