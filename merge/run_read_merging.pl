#!/usr/bin/perl -w
use strict;
use POSIX qw(tmpnam);

my ($workdir,$cycledir,$CHUNKS_FOR_MERGING,$PRIORITY) = @ARGV;


chdir($cycledir) || die "Can't go to $cycledir\n";
mkdir('CTG_MERGED') || die "Can't make CTG_MERGE directory\n";

my @this_chunk;
my @to_run;
while (my $f = <CTG_SEQ/*.ctg>) {
  push @this_chunk, $f;
  if (scalar @this_chunk >= $CHUNKS_FOR_MERGING) {
    push @to_run, [@this_chunk];
    @this_chunk=();
  }
}
if (@this_chunk) {
   push @to_run, [@this_chunk];
}


open  QSUB, ">CTG_MERGED/qsub1.sh";
print QSUB "\#!/bin/sh\n\n. /srv/sge/fedot12/common/settings.sh\n\n";

open QSUBDUMMY,">CTG_MERGED/qsub1_dummy.sh";
print QSUBDUMMY "\#!/bin/sh\n\n. /srv/sge/fedot12/common/settings.sh\n\n";

my @jobs_to_wait = ();
my $partn=0;
foreach my $chunk ( @to_run ) {
    $partn++;
    my $job_id = join( '_', "capm",$partn,get_job_id() );
    my $pfile = "capm\_$partn";
    warn "submitting $job_id\n";
    open  SH, ">CTG_MERGED/$pfile.sh";
    print SH "\#!/bin/sh\n\#\$ -S /bin/sh\n";
    print SH "mkdir /tmp/$job_id\n";
    print SH "cd /tmp/$job_id\n";
    print SH "uname -n >$pfile.host\n";
    print SH "date >>$pfile.host\n";
    print SH "cp $workdir/LIBNAMES ./\n";
    foreach my $p (@$chunk) {
      print SH "cp $cycledir/$p ./\n";
      $p=~s/.*\///;
      print SH "$workdir/normalize_coverage_and_merge_reads_multi.pl $p LIBNAMES 5 25\n";
      print SH "rm $p\n";
      print SH "mv *.reduction $cycledir/CTG_MERGED\n";
      print SH "mv *.merged $cycledir/../CTG_MERGED\n";
    }
    print SH "date >>$pfile.host\n";
    print SH "mv $pfile.host $cycledir/CTG_MERGED\n";
    print SH "cd ..\n";
    print SH "rm -R $job_id\n"; 
    close SH;
    print QSUB "qsub -P wga -o $cycledir/CTG_MERGED -e $cycledir/CTG_MERGED -N $job_id $cycledir/CTG_MERGED/$pfile.sh\n";
    push @jobs_to_wait, $job_id;
    if (scalar @jobs_to_wait > 100) {
      print QSUBDUMMY "qsub -P wga -o $cycledir/CTG_MERGED -e $cycledir/CTG_MERGED  -sync y -hold_jid ",join(",",@jobs_to_wait)," $cycledir/CTG_MERGED/dummy_cap.sh\n";
      @jobs_to_wait=();
    }
}
close QSUB;

open SH, ">CTG_MERGED/dummy_cap.sh";
print SH "\#!/bin/sh\n\#\$ -S /bin/sh\n";
close SH;

if (@jobs_to_wait) {
 print QSUBDUMMY "qsub -P wga -o $cycledir/CTG_MERGED -e $cycledir/CTG_MERGED  -sync y -hold_jid ",join(",",@jobs_to_wait)," $cycledir/CTG_MERGED/dummy_cap.sh\n";
}
close QSUBDUMMY;

system("sh $cycledir/CTG_MERGED/qsub1.sh");
system("sh $cycledir/CTG_MERGED/qsub1_dummy.sh");


sub get_job_id {
    my $id = tmpnam(); 
    $id =~ s/\/tmp\/file//;
    return $id;
}

