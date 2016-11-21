#!/usr/bin/perl -w

use strict;

my %TO_RUN;
my %JOBS;

while (my $f = <assembly1/*.sh>) {
   my ($id) = $f=~/\/(.*)\.sh/;
   next if $id=~/qsub/;
   next if $id=~/dummy/;
   my $job = `tail -n 1 $f`;
   $job=~s/\n//; $job=~s/.* //;
   $TO_RUN{$id}=$job;
   $JOBS{$job}=$id;
}

my %COMPLETED;
my %ABORTED;
my %ABORTED_ID;
my %ERR;

my $n = 0;
while (my $f = <assembly1/*.e*>) {
   my ($job) = $f=~/\/(.*)\.e/;
   next if $job=~/qsub/;
   next if $job=~/dummy/;
   $n++;
   warn "$n. $f $job\n" if $n=~/00$/;
   my $f1 = $f;
   $f1=~s/\.e/\.o/;
   my $tail = `tail $f1`;
   if ($tail=~/ALL DONE/) {
      $COMPLETED{$job}++;
   }
   else {
     $ABORTED{$job}=$tail;
     $ABORTED_ID{$job}=$f;
   }
   my $err = `cat $f`;
   if ($err=~/\S/) {
     $ERR{$job}=$err;
     $ABORTED_ID{$job}=$f;
   }
}


my %NOT_RUN;
foreach my $id (keys %TO_RUN) {
  my $job = $TO_RUN{$id};
  next if defined $COMPLETED{$job};
  next if defined $ABORTED{$job};
  $NOT_RUN{$job}++;
}

open OUT,">status.06_x1";



print OUT "Total:\t", scalar keys %TO_RUN,"\n";
print OUT "Completed:\t",scalar keys %COMPLETED,"\n";
print OUT "Aborted:\t",scalar keys %ABORTED,"\n";
print OUT "Not run:\t",scalar keys %NOT_RUN,"\n\n";
print OUT "With  errors:\t",scalar keys %ERR,"\n\n";

print OUT "Not run:\n",join("\n",keys %NOT_RUN),"\n";

print OUT "Aborted:\n";
foreach my $job (keys %ABORTED) {
  my $id = $JOBS{$job};
  print OUT "$id:\n",$ABORTED{$job},"\n\n";
}

close OUT;

open OUT,">status.06_x1.list";
foreach my $job (keys %ABORTED_ID) {
  my $id = $JOBS{$job};
  print OUT $id,"\t",$ABORTED_ID{$job},"\n";
}
foreach my $job (keys %NOT_RUN) {
  next if defined $ABORTED_ID{$job};
  my $id = $JOBS{$job};
  print OUT $id,"\t",$TO_RUN{$id},"\n";
}

close OUT;


open OUT,">status.06_err";

print OUT "Files with errors:\t", scalar keys %ERR,"\n\n";
foreach my $job (keys %ERR) {
  my $id = $JOBS{$job};
  print OUT "$id:\n",$ERR{$job},"\n\n";
}

close OUT;
