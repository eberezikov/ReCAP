#!/usr/bin/perl -w
use strict;
use POSIX qw(tmpnam);
use Config::General qw(ParseConfig);
my %CONF = ParseConfig("CONFIG");

my $RUN_HOST="/mnt/".$CONF{'RUN_HOST'};

my $PRIORITY = $ARGV[0];
   $PRIORITY = 0 unless defined $PRIORITY;
 
my @to_run;
while (my $f = <CTG_MERGED1/*.merged>) {
  push @to_run, $f;
#  last;
#  last if scalar @to_run >= 20;
}


my $workdir = `pwd`;
$workdir=~s/\n//;
$workdir=~s/\/data/$RUN_HOST/;

if (-e 'assembly1') { die "assembly1 already exists!\n\n"; }
else { mkdir 'assembly1' }

open  QSUB, ">assembly1/qsub1.sh";
print QSUB "\#!/bin/sh\n\n. /srv/sge/fedot12/common/settings.sh\n\n";

open QSUBDUMMY,">assembly1/qsub1_dummy.sh";
print QSUBDUMMY "\#!/bin/sh\n\n. /srv/sge/fedot12/common/settings.sh\n\n";

my @jobs_to_wait = ();
my $partn=0;
foreach my $p ( @to_run ) {
    $partn++;
    my $job_id = join( '_', "neb",$partn,get_job_id() );
    my $pfile = $p; $pfile=~s/.*\///; $pfile=~s/\.ctg\.merged//;
    warn "submitting $job_id\n";
    open  SH, ">assembly1/$pfile.sh";
    print SH "\#!/bin/sh\n\#\$ -S /bin/sh\n";
    print SH "mkdir /tmp/$job_id\n";
    print SH "cd /tmp/$job_id\n";
    print SH "uname -n >$pfile.host\n";
    print SH "date >>$pfile.host\n";    
    print SH "cp $workdir/$p $pfile.ctg.merged\n";
    print SH "$workdir/newbler_ctgs1.pl $pfile.ctg.merged\n";
    print SH "[ -f assembled.fa ] &&  mv assembled.fa $workdir/assembly1/$pfile.assembled.fa\n";
    print SH "[ -f assembled.fa.qual ] && mv assembled.fa.qual $workdir/assembly1/$pfile.assembled.fa.qual\n";
    print SH "[ -f assembled.reads ] && mv assembled.reads $workdir/assembly1/$pfile.assembled.reads\n";
    print SH "[ -f assembled.bl2seq ] && mv assembled.bl2seq $workdir/assembly1/$pfile.assembled.bl2seq\n";
    print SH "[ -f non-assembled.fa ] && mv non-assembled.fa $workdir/assembly1/$pfile.failed.fa\n";
    print SH "date >>$pfile.host\n";    
    print SH "mv $pfile.host $workdir/assembly1\n";
    print SH "cd ..\n";
    print SH "rm -R $job_id\n"; 
    close SH;
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

system("sh $workdir/assembly1/qsub1.sh");
system("sh $workdir/assembly1/qsub1_dummy.sh");

sub get_job_id {
    my $id = tmpnam(); 
    $id =~ s/\/tmp\/file//;
    return $id;
}

