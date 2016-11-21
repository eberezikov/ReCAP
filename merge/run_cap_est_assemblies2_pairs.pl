#!/usr/bin/perl -w
use strict;
use POSIX qw(tmpnam);

my ($workdir,$cycle,$trimming,$ASSEMBLY2_CHUNKS,$PRIORITY) = @ARGV;

my $cycledir="$workdir/MERGING/c$cycle";

$PRIORITY = 0 unless defined $PRIORITY;
   
$ASSEMBLY2_CHUNKS = 20 unless defined $ASSEMBLY2_CHUNKS;

chdir($cycledir) || die "Can't go to $cycledir\n";

if (-e 'assembly') { die "assembly already exists!\n\n"; }
else { mkdir 'assembly' }

open(F,"pairs.nofrequent.new") || die "Can't open pairs.nofrequent";
my $curr = 0;
my $n = 0;
my $p = 1;
open OUT,">assembly/p$p.grp";
while (my $l = <F>) {
  $n++;
  $curr++;
  my ($c1,$c2) = $l=~/^(\S+)\s(\S+)/;
  print OUT join("\t","s$n",$c1,$c2),"\n";
  if ($curr >= $ASSEMBLY2_CHUNKS) {
     close OUT;
     $p++;
     open OUT,">assembly/p$p.grp";
     $curr = 0;
  }
}
close OUT;

my @to_run;
while (my $f = <assembly/*.grp>) {
  push @to_run, $f;
#  last;
#  last if scalar @to_run >= 20;
}


open  QSUB, ">assembly/qsub1.sh";
print QSUB "\#!/bin/sh\n\n. /srv/sge/fedot12/common/settings.sh\n\n";

open QSUBDUMMY,">assembly/qsub1_dummy.sh";
print QSUBDUMMY "\#!/bin/sh\n\n. /srv/sge/fedot12/common/settings.sh\n\n";

my @jobs_to_wait = ();
my $partn=0;
foreach my $p ( @to_run ) {
    $partn++;
    my $job_id = join( '_', "cap2",$partn,get_job_id() );
    my $pfile = $p; $pfile=~s/.*\///; $pfile=~s/\.grp//;
    warn "submitting $job_id\n";
    open  SH, ">assembly/$pfile.sh";
    print SH "\#!/bin/sh\n\#\$ -S /bin/sh\n";
    print SH "mkdir /tmp/$job_id\n";
    print SH "cd /tmp/$job_id\n";
    print SH "uname -n >$pfile.host\n";
    print SH "date >>$pfile.host\n";    
    print SH "$workdir/cap_ctgs2_with_filtering.pl $cycledir/$p $cycledir/../CTG_MERGED $cycle $trimming\n";
    print SH "[ -f assembled.fa ] &&  mv assembled.fa $cycledir/assembly/$pfile.assembled.fa\n";
    print SH "[ -f assembled.fa.qual ] && mv assembled.fa.qual $cycledir/assembly/$pfile.assembled.fa.qual\n";
    print SH "[ -f assembled.reads ] && mv assembled.reads $cycledir/assembly/$pfile.assembled.reads\n";
    print SH "[ -f non-assembled.fa ] && mv non-assembled.fa $cycledir/assembly/$pfile.failed.fa\n";
    print SH "date >>$pfile.host\n";    
    print SH "mv $pfile.host $cycledir/assembly\n";
    print SH "cd ..\n";
    print SH "rm -R $job_id\n"; 
    close SH;
    print QSUB "qsub -P wga -o $cycledir/assembly -e $cycledir/assembly -N $job_id $cycledir/assembly/$pfile.sh\n";
    push @jobs_to_wait, $job_id;
    if (scalar @jobs_to_wait > 100) {
      print QSUBDUMMY "qsub -P wga -o $cycledir/assembly -e $cycledir/assembly  -sync y -hold_jid ",join(",",@jobs_to_wait)," $cycledir/assembly/dummy_cap.sh\n";
      @jobs_to_wait=();
    }
}
close QSUB;

open SH, ">assembly/dummy_cap.sh";
print SH "\#!/bin/sh\n\#\$ -S /bin/sh\n";
close SH;

if (@jobs_to_wait) {
 print QSUBDUMMY "qsub -P wga -o $cycledir/assembly -e $cycledir/assembly  -sync y -hold_jid ",join(",",@jobs_to_wait)," $cycledir/assembly/dummy_cap.sh\n";
}
close QSUBDUMMY;

system("sh $cycledir/assembly/qsub1.sh");
system("sh $cycledir/assembly/qsub1_dummy.sh");

sub get_job_id {
    my $id = tmpnam(); 
    $id =~ s/\/tmp\/file//;
    return $id;
}

