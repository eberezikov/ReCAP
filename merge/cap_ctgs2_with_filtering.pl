#!/usr/bin/perl -w

use strict;

my ($grp,$ctgs_dir,$cycle,$trimming) = @ARGV;

$cycle++;

my $cap3params = "-i 30 -j 31 -o 18 -s 300 -v 1 -k 0 "; #keep it same stringent at this stage
if ($trimming == 1) {
   $cap3params = "-i 30 -j 31 -o 18 -s 300 -v 1 -k 1 "; #allow trimming
}


my $totctg = `grep -c '^s' $grp`;
   $totctg=~s/\n//;

my $curr = 0;
open(GRP,$grp);


my $workdir = `pwd`;
$workdir=~s/\n//;
my ($tmpfsdir) = $workdir=~/.*\/(\S+)/;
$tmpfsdir = "/dev/shm/$tmpfsdir";
if (-e $tmpfsdir) {
  die "$tmpfsdir already exist! Aborting merge_reads\n\n";
}
else {
  mkdir $tmpfsdir;
}

chdir $tmpfsdir;

open(FINAL_FA,">$workdir/assembled.fa");
open(FINAL_QUAL,">$workdir/assembled.fa.qual");
open(FINAL_READS,">$workdir/assembled.reads");

while (my ($sname,@ctgs) = split(/\s/, <GRP>)) {
   $curr++;
   print "$curr/$totctg: $sname\n";
   my $iter = 0;
   mkdir $sname;
   chdir $sname;

   my %READS_ASSEMBLED;
   my %SIZES;
   my %PRINTED_READS;
   open REF,">$sname.ref";
   open OUT,">i$iter";
   open QUAL,">i$iter.qual";
   open CON,">i$iter.con";

   my @ctglen;
   my $readsn = 0;
   foreach my $c (@ctgs) {
     print "$curr/$totctg: $sname: reading $ctgs_dir/$c.ctg.merged\n";
     open(F,"$ctgs_dir/$c.ctg.merged") || die "Can't open $ctgs_dir/$c.ctg.merged\n";
     my $ctgid = <F>;
     my $seq = <F>;
     my $qual = <F>;
     print REF ">$c.f\n$seq";
     print OUT ">$c.f\n$seq";
     print QUAL ">$c.f\n$qual";
     $READS_ASSEMBLED{$c}=[($c)];
     push @ctglen,length($seq)-1;
     while (my $h1 = <F>) {
       last unless $h1=~/^>/;
       my $seq1 = <F>;
       my $h2 = <F>;
       my $seq2 = <F>;
       my $hq1 = <F>;
       my $qual1 = <F>;
       my $hq2 = <F>;
       my $qual2 = <F>;
       my $con = <F>;
       my $reads = <F>;
       my ($id) = $h1=~/>(\S+)\./;
       unless (defined $PRINTED_READS{$id}) {
         print OUT $h1,$seq1,$h2,$seq2;
         print QUAL $hq1,$qual1,$hq2,$qual2;
         my ($r1,$r2,$min,$max) = split(/\s/,$con);
#         if ($r1=~/^m/) { # fix sizes to old-way: wrong estimation in merge sometimes
#           $min = 80;
#           $max = 1000;
#         }
         print CON join(" ",$r1,$r2,$min,$max),"\n";
         $r1=~s/\..*//;
         $SIZES{$r1}{'min'} = $min;
         $SIZES{$r1}{'max'} = $max;
         $READS_ASSEMBLED{$r1} =[ split(/\s/, $reads) ];
         $PRINTED_READS{$id}++;
         $readsn++;
       }
     }
     close F;
   }

   close REF;
   close OUT;
   close QUAL;
   close CON;
   
   `/usr/local/blast/bin/formatdb -p F -i $sname.ref`;

   print "$curr/$totctg: $sname: assembling $readsn reads and ",scalar @ctglen," conitgs of len ",join(",", sort {$b <=> $a} @ctglen),"\n";

   while (1) {
    print "$curr/$totctg: $sname: assembly iteration $iter\n";
    `cap3 i$iter $cap3params  >i$iter.cap3`;
     my $parse_status = parse_ace($iter, $sname, \%READS_ASSEMBLED, \%SIZES);
     last if ($parse_status eq 'FINISHED');
     $iter++;
   }
   if (-e 'FAILED') {
#     `cat $ctgid.ref >>$workdir/non-assembled.fa`;
#     if ($totreadsn > 30) {
#       print "No assembly with $totreadsn reads for $ctgid!\n";
#     }
   }
   else {
      my %ASSEMBLED_SIZES;
      open(F,'final.fa');
      while (my $h = <F>) {
        my $seq = <F>;
        my ($id) = $h=~/>(\S+)/;
        $ASSEMBLED_SIZES{$id}=length($seq) - 1;
      }
      close F;
      my %ASSEMBLED_TO_SKIP;
      `/usr/local/blast/bin/blastall -p blastn -F F -d $sname.ref -i final.fa -m 9 >final.fa.bla`;
      open(BLA,"final.fa.bla");
      while (my @bla = split(/\s/,<BLA>)) {
        if ($bla[2] eq '100.00' && $bla[3] == $ASSEMBLED_SIZES{$bla[0]}) {
          $ASSEMBLED_TO_SKIP{$bla[0]}++;
        }
      }
      close BLA;

      open(F,'final.internal');
      local $/ = "\n\n";
      while (my $readblock = <F>) {
        my ($id) = $readblock=~/^(\S+)/;
        next if defined $ASSEMBLED_TO_SKIP{$id};
#        $readblock=~s/(.*\n)//; #remove header
        my $assembled_prev_ctgs = 0;
        while ($readblock=~/cap\d+\-/g) {
           $assembled_prev_ctgs++;
        }
        if ($assembled_prev_ctgs < 3) { #need at leats 3: the new id (always present) and 2 old ids
           $ASSEMBLED_TO_SKIP{$id}++;   #skip if no previously assembled contig is involved
        }
      }
      close F;
      local $/ = "\n";

      open(F,'final.fa');
      open(Q,'final.fa.qual');
      while (my $h = <F>) {
        my $seq = <F>;
        my $hq = <Q>;
        my $qual = <Q>;
        my ($id) = $h=~/>(\S+)/;
        next if defined $ASSEMBLED_TO_SKIP{$id};
        print FINAL_FA $h,$seq;
        print FINAL_QUAL $hq,$qual;
      }
      close F;
      close Q;

      open(F,'final.internal');
      local $/ = "\n\n";
      while (my $readblock = <F>) {
        my ($id) = $readblock=~/^(\S+)/;
        next if defined $ASSEMBLED_TO_SKIP{$id};
        print FINAL_READS $readblock;
      }
      close F;
      local $/ = "\n";

   }
   chdir "..";
   `rm -R $sname`;
}
close GRP;

close FINAL_FA;
close FINAL_QUAL;
close FINAL_READS;

print "\n\nALL DONE\n\n";

chdir $workdir;

`rm -R $tmpfsdir`;


##########
sub parse_ace {
  my ($iter, $origctg_id, $READS_ASSEMBLED, $SIZES) = @_;
  my $f = "i$iter";

  open ACE,"$f.cap.ace" || die "no $f.ace\n";
  local $/ = undef;
  my @ace = split(/\n\n+/,<ACE>);
  close ACE;
  local $/ ="\n";
  my $ctgname = "";
  my %CTGINFO;
  my %RDLEN;
  foreach my $block (@ace) {
   if ($block=~/^CO (\S+)/) {
     $ctgname=$1;
     $CTGINFO{$ctgname}{'CO'}=$block;
   }
   elsif ($block=~/^BQ/) {
     $CTGINFO{$ctgname}{'BQ'}=$block;
   }
   elsif ($block=~/^AF/) {
     $CTGINFO{$ctgname}{'AF'}=$block;
   }
   elsif ($block=~/^RD (\S+) (\d+)/) {
     $RDLEN{$ctgname}{$1}=$2;
   }
  }

  my %KEEP;
  my %PAIRS;
  my %REMOVE;
  my %CTGDIRECTION;
  foreach my $ctg_id (keys %CTGINFO) {
    my %DIRECTION;
    while ($CTGINFO{$ctg_id}{'AF'}=~/AF (\S+)\.(\S) (\S) (\d+)/g) {
      my ($read,$direction,$orientation,$start) = ($1,$2,$3,$4);
      $KEEP{$read}++;
      $DIRECTION{"$direction $orientation"}++;
      $PAIRS{$ctg_id}{$read}{$direction}{$orientation} = [$start, $start+$RDLEN{$ctg_id}{"$read\.$direction"}];
    }
    my $plus_direction = 0;
    $plus_direction+=$DIRECTION{"f U"} if defined $DIRECTION{"f U"};
    $plus_direction+=$DIRECTION{"r C"} if defined $DIRECTION{"r C"};
    my $minus_direction = 0;
    $minus_direction+=$DIRECTION{"f C"} if defined $DIRECTION{"f C"};
    $minus_direction+=$DIRECTION{"r U"} if defined $DIRECTION{"r U"};

    if ($minus_direction > $plus_direction) {
      $CTGDIRECTION{$ctg_id}= -1;
    }
    else {
      $CTGDIRECTION{$ctg_id}= 1;
    }
  }
 
  open(CON,">$f.con.myresults");
  foreach my $ctg_id (keys %PAIRS) {
   foreach my $id (keys %{$PAIRS{$ctg_id}}) {
    if ($id=~/^cap/) {
      print CON "$id\toriginal ctg\tkeep\n";
      next;
    }
    unless (defined $PAIRS{$ctg_id}{$id}{'f'} && defined $PAIRS{$ctg_id}{$id}{'r'}) {
      $REMOVE{$id}++;
      delete $KEEP{$id};
      print CON "$id\tno pair\tremoved\n";
      next;
    }

#   only checking orientation and distance if both in the same contig
     my @dist;
     if ($CTGDIRECTION{$ctg_id} == 1) {
      if (defined $PAIRS{$ctg_id}{$id}{'f'}{'C'} || defined $PAIRS{$ctg_id}{$id}{'r'}{'U'}) {
          $REMOVE{$id}++;
          delete $KEEP{$id};
          print CON "$id\twrong strand\tremoved\n";
          next;
      }
      else {
           push @dist, @{$PAIRS{$ctg_id}{$id}{'f'}{'U'}}, @{$PAIRS{$ctg_id}{$id}{'r'}{'C'}} if (defined $PAIRS{$ctg_id}{$id}{'f'}{'U'} && defined $PAIRS{$ctg_id}{$id}{'r'}{'C'});
      }
     } 
     else {
      if (defined $PAIRS{$ctg_id}{$id}{'f'}{'U'} || defined $PAIRS{$ctg_id}{$id}{'r'}{'C'}) {
          $REMOVE{$id}++;
          delete $KEEP{$id};
          print CON "$id\twrong strand\tremoved\n";
          next;
      }
      else {
          push @dist, @{$PAIRS{$ctg_id}{$id}{'f'}{'C'}}, @{$PAIRS{$ctg_id}{$id}{'r'}{'U'}} if (defined $PAIRS{$ctg_id}{$id}{'f'}{'C'} && defined $PAIRS{$ctg_id}{$id}{'r'}{'U'});
      }
    }
    if (@dist) {
     @dist = sort {$a <=> $b} @dist;
     my $len = $dist[3] - $dist[0] + 1;
     if ($len < $SIZES->{$id}{'min'} || $len > $SIZES->{$id}{'max'}) {
        $REMOVE{$id}++;
        delete $KEEP{$id};
        print CON join("\t",$id,$len,$SIZES->{$id}{'min'},$SIZES->{$id}{'max'},'removed wrong size'),"\n";
     }
     else {
        print CON join("\t",$id,$len,$SIZES->{$id}{'min'},$SIZES->{$id}{'max'},'satisfied'),"\n";
     }
    }
   }
  }
  close CON;

  if (scalar keys %KEEP == 0) {
     print "Failed assembly, no reads remained\n";
     `touch FAILED`;
     return 'FINISHED';
  }

  if (scalar keys %REMOVE > 0) {
    $iter++;
    my $fnext="i$iter";
    open F,"$f";
    open QUAL,"$f.qual";
    open OUT,">$fnext";
    open OUTQUAL,">$fnext.qual";
    while (my $h = <F>) {
      my $seq = <F>;
      my $hq = <QUAL>;
      my $qual = <QUAL>;
      my ($id) = $h=~/>(\S+)\./;
      next unless defined $KEEP{$id};
      print OUT $h,$seq;
      print OUTQUAL $hq,$qual;
    }
    close F;
    close QUAL;
    close OUT;
    close OUTQUAL;
   
    open F,"$f.con";
    open OUT,">$fnext.con";
    while (my $l = <F>) {
      my ($id) = $l=~/^(\S+)\./;
      next unless defined $KEEP{$id};
      print OUT $l;
    }
    close F;
    close OUT;
    return "Do another iteration";
  }
  else {
    open OUT,">final.fa";
    open OUTQUAL,">final.fa.qual";
    open INTERNALLIST,">final.internal";
    my $this_ctg_n = 0;
    my ($orign) = $origctg_id=~/s(\d+)/;
    foreach my $ctg_id (keys %CTGINFO) {
      $this_ctg_n++;
      my $this_ctg_id = 'cap'.$cycle.'-'.$orign.'-'.$this_ctg_n;
      my $ctgseq = $CTGINFO{$ctg_id}{'CO'};
      my $ctgqual = $CTGINFO{$ctg_id}{'BQ'};
      $ctgseq=~s/CO.*\n//;
      $ctgseq=~s/\n//g;
      $ctgseq=~s/\*//g;
      $ctgqual=~s/BQ.*\n//;
      $ctgqual=~s/\n/ /g; $ctgqual=~s/\s+/ /g; $ctgqual=~s/\s+$//;
      my $ctgseq_len = length($ctgseq);
      my %READS_TO_PRINT;
      foreach my $id (keys %{$PAIRS{$ctg_id}}) {
        next unless defined $KEEP{$id};
        foreach my $rid (@{$READS_ASSEMBLED->{$id}}) {
          $READS_TO_PRINT{$rid}++;
        }
      }
      my $nreads = scalar keys %READS_TO_PRINT;
      print INTERNALLIST join(" ",$this_ctg_id,$ctgseq_len,$nreads,sprintf("%.2f",$nreads/$ctgseq_len)),"\n";
      print INTERNALLIST join("\n",keys %READS_TO_PRINT),"\n\n";
      if ($CTGDIRECTION{$ctg_id} == -1) {
         $ctgseq = revcomp($ctgseq);
         $ctgqual = join(" ",reverse(split(/\s+/,$ctgqual)));
      }
      print OUT ">$this_ctg_id ",length($ctgseq),"\n$ctgseq\n";
      print OUTQUAL ">$this_ctg_id\n$ctgqual\n";
    }
    close OUT;
    close OUTQUAL;
    close INTERNALLIST;
    print "All done - finishing at this iteration\n";
    return 'FINISHED';
  }
}

##########


##########
sub revcomp {
  my $str = shift;
  $str=reverse($str);
  $str=~tr/ATGCatgc/TACGtacg/;
  return $str;
}
##########
