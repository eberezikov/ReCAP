#!/usr/bin/perl -w

use strict;


my $MINREADS = 1;
my $MIN_FRACTION_FOR_CLUSTERING = 0.50;


my %CLUSTERS;
my %ASSIGNED;
my %CLUSTERREADS;

my $n = 0;
open(F,"FINAL/contributions_merged.txt");
while (my ($id1,$totreads,$selfreads,$ctgn,@contributors) = split(/\s/,<F>)) {
    $n++;
#    last if $n > 10;
#    warn "Considering $id1 $totreads $selfreads\n";
#    warn join("\t","",$id1,$totreads,$selfreads,$ctgn,@contributors),"\n";
    my $reads_for_clustering = $selfreads*$MIN_FRACTION_FOR_CLUSTERING;
    my $cluster_id;
    foreach my $c (@contributors) {
       my ($id2,$r) = split(/\:/,$c);
       next if  $r < $reads_for_clustering;
       if (defined $ASSIGNED{$id2}) {
         $cluster_id = $ASSIGNED{$id2};
 #       warn "assigning to the existing cluster $cluster_id\n";
         last;
       }
    }
    unless (defined $cluster_id) {
      $cluster_id = $id1;
      $ASSIGNED{$id1}=$cluster_id;
#      warn "starting new cluster with $cluster_id\n";
    }
    foreach my $c (@contributors) {
       my ($id2,$r) = split(/\:/,$c);
#       warn "  considering $id2,$r\n";
       next if  $r < $reads_for_clustering;
       unless (defined $CLUSTERS{$cluster_id}{$id2}) {  #to skip repeated contrib counts
         $CLUSTERS{$cluster_id}{$id2}=$r ;
         $CLUSTERREADS{$cluster_id}+=$r;
       }
       $ASSIGNED{$id2}=$cluster_id;
#       warn "   added to cluster $cluster_id: $id2\n";
#       warn "   $id2 assigned to cluster $cluster_id\n";
    }
}

open OUT1,">FINAL/clusters.txt";
open OUT2,">FINAL/clustered_transcripts.txt";
foreach my $cl (sort {$CLUSTERREADS{$b} <=> $CLUSTERREADS{$a}} keys %CLUSTERREADS) {
    print OUT1 join("\t",$cl,$CLUSTERREADS{$cl}, scalar keys %{$CLUSTERS{$cl}});
    foreach my $tr (sort {$CLUSTERS{$cl}{$b} <=> $CLUSTERS{$cl}{$a}} keys %{$CLUSTERS{$cl}}) {
       print OUT1 "\t","$tr:$CLUSTERS{$cl}{$tr}";
       print OUT2 $tr,"\t",$cl,"\n";
    }
    print OUT1 "\n";
}
close OUT1;
close OUT2;