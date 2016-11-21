#!/usr/bin/perl -w

use strict;

my ($ctgs) = @ARGV;

#my $cap3params = "-i 30 -j 31 -o 18 -s 300 -v 1 -k 0 ";

my $totctg = `grep -c contig $ctgs`;
   $totctg=~s/\n//;

my $curr = 0;
open(BLOCK,$ctgs);
local $/ = "\n\n";


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

while (my $block = <BLOCK>) {
   $curr++;
   $block=~s/\n\n//;
   my ($ctgid,$ctgseq,$ctgqual,@reads) = split(/\n/,$block);
   my $ctglen = length($ctgseq);
   print "$curr/$totctg: $ctgid ($ctglen)\n";
   my $iter = 0;
   mkdir $ctgid;
   chdir $ctgid;
   open REF,">$ctgid.ref";
   print REF ">$ctgid\n$ctgseq\n";
   close REF;
   `/usr/local/blast/bin/formatdb -p F -i $ctgid.ref`;
   my $totreadsn = (scalar @reads)/10;
   print "$curr/$totctg: $ctgid ($ctglen): reads = $totreadsn\n";
   open OUT,">i$iter";
   open QUAL,">i$iter.qual";
   open CON,">i$iter.con";
   my %READS_ASSEMBLED;
   my %SIZES;
   my $last_read = $#reads - 8;
   for (my $i=0; $i < $last_read; $i+=10) {
      my ($id1,$direction1) = $reads[$i]=~/>(\S+)\.(\S+)/;
      my ($lib1) = $id1=~/^(L\d+|m)/;
      $direction1=uc($direction1);
      my $newbler_info1 = " template=$id1 dir=$direction1 library=$lib1";
      $reads[$i].= $newbler_info1;
      $reads[$i+4].= $newbler_info1;
      my ($id2,$direction2) = $reads[$i+2]=~/>(\S+)\.(\S+)/;
      my ($lib2) = $id2=~/^(L\d+|m)/;
      $direction2=uc($direction2);
      my $newbler_info2 = " template=$id2 dir=$direction2 library=$lib2";
      $reads[$i+2].= $newbler_info2;
      $reads[$i+6].= $newbler_info2;
      print OUT join("\n",@reads[$i..$i+3]),"\n";
      print QUAL join("\n",@reads[$i+4..$i+7]),"\n";
      my ($r1,$r2,$min,$max) = split(/\s/,$reads[$i+8]);
#      if ($r1=~/^m/) { # fix sizes to old-way: wrong estimation in merge sometimes
#        $min = 80;
#        $max = 1000;
#      }
      print CON join(" ",$r1,$r2,$min,$max),"\n";
      $r1=~s/\..*//;
      $SIZES{$r1}{'min'} = $min;
      $SIZES{$r1}{'max'} = $max;
      $READS_ASSEMBLED{$r1} =[ split(/\s/, $reads[$i+9]) ];
   }
   while (1) {
    print "$curr/$totctg: $ctgid ($ctglen): assembly iteration $iter\n";
    `/opt/454/bin/runAssembly -cdna -urt -nobig -ace -o i$iter.nb i$iter >i$iter.log 2>&1`;
     my $parse_status = parse_ace($iter, $ctgid, \%READS_ASSEMBLED, \%SIZES);
     last if ($parse_status eq 'FINISHED');
     $iter++;
   }
   if (-e 'FAILED') {
     `cat $ctgid.ref >>$workdir/non-assembled.fa`;
     if ($totreadsn > 30) {
       print "No assembly with $totreadsn reads for $ctgid!\n";
#       chdir $workdir;
#       `rm -R $tmpfsdir`;
#       warn "Terminated due to suspicious lack of assembly for $ctgid ($totreadsn reads).\n";
#       exit();
     }
   }
   else {
     `cat final.fa >>$workdir/assembled.fa`;
     `cat final.fa.qual >>$workdir/assembled.fa.qual`;
     `cat final.internal >>$workdir/assembled.reads`;
     `cat bl2seq.hit >>$workdir/assembled.bl2seq`;
   }
   chdir "..";
   `rm -R $ctgid`;
}
close BLOCK;

print "\n\nALL DONE\n\n";

chdir $workdir;

`rm -R $tmpfsdir`;


##########
sub parse_ace {
  my ($iter, $origctg_id, $READS_ASSEMBLED, $SIZES) = @_;
  my $f = "i$iter";

  unless (-e "i$iter.nb/454Isotigs.ace") {
     print "Failed assembly, no ace file generated\n";
     `touch FAILED`;
     return 'FINISHED';
  }
  open(ACE,"i$iter.nb/454Isotigs.ace") || die "no i$iter.nb/454Isotigs.ace\n";
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
   elsif ($block=~/^RD (\S+?\.\S)\S* (\d+)/) {
     $RDLEN{$ctgname}{$1}=$2;
   }
  }

  my %KEEP;
  my %PAIRS;
  my %REMOVE;
  my %CTGDIRECTION;
  foreach my $ctg_id (keys %CTGINFO) {
    my %DIRECTION;
    while ($CTGINFO{$ctg_id}{'AF'}=~/AF (\S+?)\.(\S)\S* (\S) (\d+)/g) {
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
    my ($orign) = $origctg_id=~/-(\d+)/;
    foreach my $ctg_id (keys %CTGINFO) {
      $this_ctg_n++;
      my $this_ctg_id = 'cap1-'.$orign.'-'.$this_ctg_n;
      my $ctgseq = $CTGINFO{$ctg_id}{'CO'};
      my $ctgqual = $CTGINFO{$ctg_id}{'BQ'};
      $ctgseq=~s/CO.*\n//;
      $ctgseq=~s/\n//g;
      $ctgseq=~s/\*//g;
      $ctgqual=~s/BQ.*\n//;
      $ctgqual=~s/\n/ /g; $ctgqual=~s/\s+/ /g; $ctgqual=~s/\s+$//;
      my $ctgseq_len = length($ctgseq);
      my %READS_TO_PRINT;
      foreach my $id (keys %KEEP) {
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
    `/usr/local/blast/bin/blastall -p blastn -F F -d $origctg_id.ref -i final.fa -m 8 -W 30 >bl2seq.hit`;
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
