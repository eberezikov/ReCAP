#!/usr/bin/perl -w
use strict;
use POSIX qw(tmpnam);

my ($readsdir, $workdir, $K, $P, $NM, $PRIORITY) = @ARGV;

$PRIORITY = 0 unless defined $PRIORITY;

my $GENOME = "$workdir/FINAL/final_all.fa";
my $outdir = "$workdir/FINAL/map";

my @to_run;
while (my $f = <$readsdir/*/*.R1.fq>) {
  push @to_run, $f;
#  last;
}

open  QSUB, ">$workdir/FINAL/bwt.sh";
print QSUB "\#!/bin/sh\n\n. /srv/sge/fedot12/common/settings.sh\n\n";


my @jobs_to_wait = ();
my $partn=0;
foreach my $p ( @to_run ) {
    $partn++;
    my $job_id = join( '_', "bwt",$partn,get_job_id() );
    my $pfile = $p; $pfile=~s/.*\///; $pfile=~s/\.R1\.fq//;
    my ($pdir) = $p=~/$readsdir\/(\S+)\//;
    mkdir "$outdir/$pdir" unless -e "$outdir/$pdir";
    warn "submitting $pfile\n";
    open  SH, ">$outdir/$pdir/$job_id.sh";
    print SH "\#!/bin/sh\n\#\$ -S /bin/sh\n";
    print SH "mkdir /tmp/$job_id\n";
    print SH "cd /tmp/$job_id\n";
    print SH "uname -n >$pfile.host\n";
    print SH "cp $readsdir/$pdir/$pfile.R1.fq $readsdir/$pdir/$pfile.R2.fq ./\n";
    print SH "date >>$pfile.host\n";
    print SH "$workdir/run_bowtie_global_filter_best_and_distance_and_strand.pl $GENOME $pfile.R1.fq $pfile.R1.sam $pfile.R1.err $P $K $NM R\n";
    print SH "$workdir/run_bowtie_global_filter_best_and_distance_and_strand.pl $GENOME $pfile.R2.fq $pfile.R2.sam $pfile.R2.err $P $K $NM F\n";
    print SH "$workdir/filter_nonmapped_with_size_check.pl $workdir/LIBNAMES $GENOME $NM $K\n";
    print SH "date >>$pfile.host\n";
    print SH "cp *.used *.nonmapped.fq $pfile.host $pfile.*.err $outdir/$pdir\n";
    print SH "cp pairs $outdir/$pdir/$pfile.map\n";
    print SH "cp links $outdir/$pdir/$pfile.lnk\n";
    print SH "cp freq $outdir/$pdir/$pfile.frq\n";
    print SH "cd ..\n";
    print SH "rm -R /tmp/$job_id\n"; 
    close SH;
    print QSUB "qsub -P wga -pe threaded $P -o $outdir/$pdir -e $outdir/$pdir -N $job_id $outdir/$pdir/$job_id.sh\n";
    push @jobs_to_wait, $job_id;
}

open SH, ">$workdir/FINAL/dummy_bwt.sh";
print SH "\#!/bin/sh\n\#\$ -S /bin/sh\n";
close SH;

print QSUB "qsub -P wga -o $workdir/FINAL -e $workdir/FINAL  -sync y -hold_jid ",join(",",@jobs_to_wait)," $workdir/FINAL/dummy_bwt.sh";
close QSUB;

system("sh $workdir/FINAL/bwt.sh");

sub get_job_id {
    my $id = tmpnam(); 
    $id =~ s/\/tmp\/file//;
    return $id;
}
