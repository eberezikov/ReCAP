#!/usr/bin/perl -w

use strict;
use Parallel::ForkManager;

my ($CTG_SEQ,$RAW_READS) = @ARGV;

my $pm = Parallel::ForkManager->new(21);

warn "Reading ctg length\n";
my %LEN;
open(F,"FINAL/final_all.fa");
while (my $h = <F>) {
  my $seq = <F>;
  my ($id,$len) = $h=~/>(\S+) (\d+)/;
  $LEN{$id}=$len;
}
close F;

my @libs = <FINAL/map/*>;

warn "Will process these libs: ",join(" ",@libs),"\n";

foreach my $libdir (@libs) {
 $pm->start and next; 
  my %CTGS;
  my @reads;
  open LOG,">$libdir/priority.log";
  { my $ofh = select LOG; #for immediate write
    $| = 1;
    select $ofh;
  }
  while (my $f = <$libdir/*.mapped_links>) {
    my $tot = `wc -l $f`;
    $tot=~s/\s.*\n//;
    my $curr = 0;
    open(F,$f);
    while (my ($r,@ctgs) = split(/\s/,<F>)) {
      $curr++;
#      last if $curr > 10000;
      if ($curr=~/00000$/) {
        print LOG "$f. $curr/$tot (",sprintf("%.2f",$curr/$tot*100),"%)\n";
      }
      next unless @ctgs;
      $r=~s/\/.*//;
      my ($lib,$n) = $r=~/L(\d+)r(\d+)/;
      foreach my $c (@ctgs) {
        $c=~s/\/.*//;
        push @{$CTGS{$c}},$r;
        $reads[$lib][$n].="$c ";
      }
    }
    close F;
  }

  print LOG "Reading internal reads\n";
  while (my $f = <$libdir/*.mapped_reads>) {
    open(F,$f);
    my $tot = `wc -l $f`;
    $tot=~s/\s.*\n//;
    my $curr = 0;
    while (my @l = split(/\s/,<F>)) {
      push @{$CTGS{$l[2]}},$l[1];
      my ($lib,$n) = $l[1]=~/L(\d+)r(\d+)/;
      $reads[$lib][$n].="$l[2] ";
      $curr++;
#      last if $curr > 10000;
      if ($curr=~/00000$/) {
        print LOG "$f. $curr/$tot (",sprintf("%.2f",$curr/$tot*100),"%)\n";
      }
    }
    close F;
  }

  print LOG "Prioritizing contigs\n";

  my %ASSIGNED_CTGS;
  my $totassigned = 0;
  my $totctgs = scalar keys %CTGS;

  my %FREQ;
  my %PREVIOUS_FREQ;
  my %FREQ_DISTRIBUTION;
  my %FREQ_CONTRIBUTION;
  my @ranges = (10000,8000,5000,3000,1000,900,800,700,600,500,400,300,200,100,95,90,85,80,95,70,65,60,55,50,45,40,38,36,34,32,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0);
  my %RANGE;
  foreach my $c (keys %CTGS) {
    $FREQ{$c} = scalar @{$CTGS{$c}};
    $PREVIOUS_FREQ{$c} = $FREQ{$c};
    my $range = 0;
    foreach my $r (@ranges) {
      if ($FREQ{$c} > $r) {
         $range = $r;
         $RANGE{$r}{$c} = $FREQ{$c};
         last;
      }
    }
    $FREQ_DISTRIBUTION{$range}++;
  }
  
  for my $range (sort {$b <=> $a} keys %FREQ_DISTRIBUTION) {
     print LOG "Frequency > $range : $FREQ_DISTRIBUTION{$range}\n";
  }
  
  my %ORDERED;
  print LOG "Sorting in ranges\n";
  foreach my $r (@ranges) {
     print LOG "Sorting in range $r\n";
     @{$ORDERED{$r}} = sort { $RANGE{$r}{$b} <=> $RANGE{$r}{$a} || $LEN{$b} <=> $LEN{$a} } keys %{$RANGE{$r}};
  }

  
  my $current_range = $ranges[0];
  my @ordered = @{$ORDERED{$current_range}};
  @{$ORDERED{$current_range}}=();
  my @lower_range;
#  my @reverse_ranges = reverse(@ranges);

  while ($totassigned < $totctgs) {

    my $resort = 0;
    if (scalar @ordered == 0 ) {
      print LOG "No more contigs in range $current_range. Recursing\n";
      foreach my $r (@ranges) {
        if (@{$ORDERED{$r}}) {
           @ordered = @{$ORDERED{$r}};
           print LOG "ADDED ",scalar keys @{$ORDERED{$r}}," contigs of range $r. Current range is $r.\n";
           @{$ORDERED{$r}} = ();
           $current_range = $r;
           $resort=1;
           last;
        }
      }
     }

     if (scalar @ordered == 0) {
       print LOG "No more contigs to process. Finishing.\n";
       last;
     }

     print LOG "$totassigned/$totctgs. Finding the most frequent ctg\n";
     my $prev_top = $ordered[0];
     $PREVIOUS_FREQ{$prev_top}=0 unless defined $PREVIOUS_FREQ{$prev_top};
     if ($FREQ{$prev_top} != $PREVIOUS_FREQ{$prev_top}) { #need to resort, since the potential higest was affected
       print LOG "  resorting frequencies for ",scalar @ordered," contigs. Current range is $current_range.\n";
       my @ordered_tmp;
       if ($current_range > 0) {
          @ordered_tmp = sort { $FREQ{$b} <=> $FREQ{$a} || $LEN{$b} <=> $LEN{$a} } @ordered;
       }
       else {
          @ordered_tmp = sort { $LEN{$b} <=> $LEN{$a} } @ordered;       
       }
       @ordered=();
       foreach my $id (@ordered_tmp) {
         if ($FREQ{$id} > $current_range) {
            push @ordered,$id;
          }
          else {
            foreach my $r (@ranges) {
              if ($FREQ{$id} > $r) {
                 push @{$ORDERED{$r}}, $id;
                 last;
              }
            }
          }
        }
        if (defined $ordered[1]) {
          $PREVIOUS_FREQ{$ordered[1]}=$FREQ{$ordered[1]}; #updating for next top for the next round
        }
     }

    if (scalar @ordered == 0 ) {
      print LOG "No more contigs in range $current_range. Recursing\n";
      foreach my $r (@ranges) {
        if (@{$ORDERED{$r}}) {
           @ordered = @{$ORDERED{$r}};
           print LOG "ADDED ",scalar keys @{$ORDERED{$r}}," contigs of range $r. Current range is $r.\n";
           @{$ORDERED{$r}} = ();
           $current_range = $r;
           $resort=1;
           last;
        }
      }
     }

     if (scalar @ordered == 0) {
       print LOG "No more contigs to process. Finishing.\n";
       last;
     }

     my $top = $ordered[0];

     foreach my $r (@ranges) {
       if ($FREQ{$top} <= $r && @{$ORDERED{$r}}) {
          push @ordered, @{$ORDERED{$r}};
          print LOG "ADDED ",scalar keys @{$ORDERED{$r}}," contigs of range $r. Current range is $r.\n";
          @{$ORDERED{$r}} = ();
          $current_range = $r;
          $resort=1;
          last;
       }
     }

     if ($resort == 1) { #sort again if more were added
        print LOG "  resorting frequencies for ",scalar @ordered," contigs. Current range is $current_range.\n";
        my @ordered_tmp;
        if ($current_range > 0) {
           @ordered_tmp = sort { $FREQ{$b} <=> $FREQ{$a} || $LEN{$b} <=> $LEN{$a} } @ordered;
        }
        else {
          @ordered_tmp = sort { $LEN{$b} <=> $LEN{$a} } @ordered;
        }
        @ordered=();
        foreach my $id (@ordered_tmp) {
          if ($FREQ{$id} > $current_range) {
             push @ordered,$id;
           }
           else {
             foreach my $r (@ranges) {
               if ($FREQ{$id} > $r) {
                  push @{$ORDERED{$r}}, $id;
                  last;
               }
             }
           }
        }
        if (defined $ordered[1]) {
           $PREVIOUS_FREQ{$ordered[1]}=$FREQ{$ordered[1]}; #updating for next top for the next round
        }
     }

     my $ctg_id = shift @ordered;
     my $ctg_freq = $FREQ{$ctg_id};
     print LOG "  $ctg_id -> $ctg_freq reads\n";
     $ASSIGNED_CTGS{$ctg_id}=$ctg_freq;
     my $curr = 0;
     foreach my $r (@{$CTGS{$ctg_id}}) {
        $curr++;
        if ($curr=~/0000$/) { 
           print LOG "$totassigned/$totctgs. Processed $curr/$ctg_freq reads (",sprintf("%.2f",$curr/$ctg_freq*100),"%)\n";
        }
        my ($lib,$n) = $r=~/L(\d+)r(\d+)/;
        next if $reads[$lib][$n]=~/^1/;
        foreach my $c (split(/\s+/,$reads[$lib][$n])) {
          next unless defined $FREQ{$c};
          $FREQ{$c}--;
          $FREQ_CONTRIBUTION{$c}{$ctg_id}++;
        }
        $reads[$lib][$n]="1 $ctg_id";
     }

     $totassigned = scalar keys %ASSIGNED_CTGS;

  }

  print LOG "Writing prioritized contigs\n";
  open OUT,">$libdir/priority.txt";
  foreach my $id (sort {$ASSIGNED_CTGS{$b} <=> $ASSIGNED_CTGS{$a}} keys %ASSIGNED_CTGS) {
     print OUT join("\t",$id,$LEN{$id},$ASSIGNED_CTGS{$id}),"\n";
  }
  close OUT;
  
  print LOG "Writing ctg contributions\n";
  open OUT,">$libdir/contributions.txt";
  foreach my $id (sort keys %FREQ_CONTRIBUTION) {
    print OUT $id;
    foreach my $c (sort {$FREQ_CONTRIBUTION{$id}{$b} <=> $FREQ_CONTRIBUTION{$id}{$a}} keys %{$FREQ_CONTRIBUTION{$id}}) {
       print OUT "\t",$c,':',$FREQ_CONTRIBUTION{$id}{$c};
    } 
    print OUT "\n";
  }

  open OUT,">$libdir/noreads_ctgs.txt";
  foreach my $id (sort keys %LEN) {
     next if defined $ASSIGNED_CTGS{$id};
     print OUT $id,"\n";
  }
  close OUT;
  close LOG;
 $pm->finish;
} 

$pm->wait_all_children();

