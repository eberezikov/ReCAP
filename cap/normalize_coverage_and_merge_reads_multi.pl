#!/usr/bin/perl -w

use strict;

my ($ctgs,$libnames,$max_reads_per_position,$window_for_merging) = @ARGV;

my $cap3params = "-i 30 -j 31 -o 18 -s 300 -v 1 -k 0";

$max_reads_per_position = 3 unless defined $max_reads_per_position;
$window_for_merging = 30 unless defined $window_for_merging;

my $max_reads_per_batch = 5000000;
$max_reads_per_batch = $max_reads_per_batch*8; #fatsq pe


my %SIZES;
my %ORIENTATION;
open(F,$libnames);
while (my @l = split(/\s/,<F>)) {
   $SIZES{$l[1]}{'min'}=$l[2];
   $SIZES{$l[1]}{'max'}=$l[3];
   $ORIENTATION{$l[1]}=$l[4];
}
close F;


my $ctgs_to_process = `grep -c \'^\$\' $ctgs`;
$ctgs_to_process=~s/\n//;
$ctgs_to_process++;

my $current = 0;

open(BLOCK,$ctgs);

open OUT,">$ctgs.merged";
open STAT,">$ctgs.reduction";

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


while (my $ctgid = <BLOCK>) {
  my $ctgseq = <BLOCK>;
  my $ctgqual = <BLOCK>;
  print OUT $ctgid,$ctgseq,$ctgqual;
  my $ctglen = length($ctgseq)-1;
     $ctgid=~s/\n//;
  my ($unqid) = $ctgid=~/^\S+?-(\S+)/;
  $current++;
  print "$current/$ctgs_to_process: $ctgid ($ctglen): making reference\n";
  open REF,">$ctgid.ref";
  print REF ">$ctgid\n$ctgseq";
  close REF;
  `/usr/local/bowtie/2.2.4/bowtie2-build $ctgid.ref $ctgid.ref 2>/dev/null`;
  my @fastq_batches;
  my $fastq_part = 1;
  print "$current/$ctgs_to_process: $ctgid ($ctglen): making fastq batch $fastq_part\n";
  open R,">in.$fastq_part.fastq";
  my $totreadsn = 0;
  my $curreadsn = 0;
  my $pairn = 0;
  while (my $l = <BLOCK>) {
      last if $l=~/^\n/;
      print R $l;
      $totreadsn++;
      $curreadsn++;
      if ($curreadsn >= $max_reads_per_batch) {
        close R;
        push @fastq_batches, "in.$fastq_part.fastq";
        $fastq_part++;
        open R,">in.$fastq_part.fastq";
        print "$current/$ctgs_to_process: $ctgid ($ctglen): making fastq batch $fastq_part\n";
        $curreadsn=0;
      }
  }

  close R;

  if (-s "in.$fastq_part.fastq" > 0) {
     push @fastq_batches, "in.$fastq_part.fastq";
  }
  $totreadsn = $totreadsn/8;
  my $totbatches = scalar @fastq_batches;
  print "$current/$ctgs_to_process: $ctgid ($ctglen): reads = $totreadsn, batches = $totbatches\n";
  my $currbatch = 0;

  my %ALL_READS_PER_POS1;
  my %ALL_READS_PER_POS2;

  foreach my $fastq (@fastq_batches) {
    $currbatch++;
    print "$current/$ctgs_to_process: $ctgid ($ctglen), batch $currbatch/$totbatches: reading fastq\n";
    my %INPUTFASTQ;
    open(FQ,$fastq);
    while (my $h1 = <FQ>) {
      my $seq1 = <FQ>;
      my $hq1 = <FQ>;
      my $qual1 = <FQ>;
      my $h2 = <FQ>;
      my $seq2 = <FQ>;
      my $hq2 = <FQ>;
      my $qual2 = <FQ>;
      my ($fqid) = $h1=~/\@(\S+)\//;
      $INPUTFASTQ{$fqid} = $h1.$seq1.$hq1.$qual1.$h2.$seq2.$hq2.$qual2;
    }
    close(FQ);
    my %POS1;
    my %POS2;
    print "$current/$ctgs_to_process: $ctgid ($ctglen), batch $currbatch/$totbatches: mapping reads\n";
    open (BWT,"/usr/local/bowtie/2.2.4/bowtie2 -x $ctgid.ref -q -U $fastq --end-to-end  --reorder --sam-nohead --sam-nosq 2>bwt.err |");
    while (my @l1 = split(/\t/,<BWT>)) {
      my @l2 = split(/\t/,<BWT>);
      my ($id) = $l1[0]=~/(\S+)\//;
      my ($p1,$p2) = ($l1[3],$l2[3]);
      unless ($p1 == 0) { push @{$POS1{$p1}}, $id; }
      unless ($p2 == 0) { push @{$POS2{$p2+length($l2[9])}}, $id; }
    }
    print "$current/$ctgs_to_process: $ctgid ($ctglen), batch $currbatch/$totbatches: parsing ",scalar keys %POS1," positions\n";
    foreach my $p (sort {$a <=> $b} keys %POS1) {
      my (%FASTQ,%QUALSUM);
      my $nid = 0;
      foreach my $id (@{$POS1{$p}}) {
          $nid++;
          last if  $nid > 500; #consider max 500 reads per position per batch for quality calculation, to speed up calculations
          my @fq = split(/\n/,$INPUTFASTQ{$id});
          my $qual1_sum = count_qual($fq[3]);
          my $qual2_sum = count_qual($fq[7]);
          $QUALSUM{$id}=$qual1_sum+$qual2_sum;
      }
      my @ids = sort {$QUALSUM{$b} <=> $QUALSUM{$a}} keys %QUALSUM;
      for my $i (0..$max_reads_per_position-1) {
         last unless defined $ids[$i];
         push @{$ALL_READS_PER_POS1{$p}},[$QUALSUM{$ids[$i]},$INPUTFASTQ{$ids[$i]}];
      }
    }
    foreach my $p (sort {$b <=> $a} keys %POS2) {
      my (%FASTQ,%QUALSUM);
      my $nid = 0;
      foreach my $id (@{$POS2{$p}}) {
          $nid++;
          last if  $nid > 500; #consider max 500 reads per position per batch for quality calculation, to speed up calculations
          my @fq = split(/\n/,$INPUTFASTQ{$id});
          my $qual1_sum = count_qual($fq[3]);
          my $qual2_sum = count_qual($fq[7]);
          $QUALSUM{$id}=$qual1_sum+$qual2_sum;
      }
      my @ids = sort {$QUALSUM{$b} <=> $QUALSUM{$a}} keys %QUALSUM;
      for my $i (0..$max_reads_per_position-1) {
         last unless defined $ids[$i];
         push @{$ALL_READS_PER_POS2{$p}},[$QUALSUM{$ids[$i]},$INPUTFASTQ{$ids[$i]}];
      }
    }
  }


  #resorting to merge batches
  my %BEST_READS_PER_POS1;
  my %BEST_READS_PER_POS2;
  foreach my $p (keys %ALL_READS_PER_POS1) {
     my @fq = sort {$b->[0] <=> $a->[0]} @{$ALL_READS_PER_POS1{$p}};
     for my $i (0..$max_reads_per_position-1) {
         last unless defined $fq[$i];
         push @{$BEST_READS_PER_POS1{$p}},$fq[$i]->[1];
     }
  }
  foreach my $p (keys %ALL_READS_PER_POS2) {
     my @fq = sort {$b->[0] <=> $a->[0]} @{$ALL_READS_PER_POS2{$p}};
     for my $i (0..$max_reads_per_position-1) {
         last unless defined $fq[$i];
         push @{$BEST_READS_PER_POS2{$p}},$fq[$i]->[1];
     }
  }

  #aggregating max best per position per fixed window
  my %POS1;
  foreach my $p (keys %BEST_READS_PER_POS1) {
     my $p1 = int($p/$window_for_merging)+1;
     push @{$POS1{$p1}},@{$BEST_READS_PER_POS1{$p}};
  }
  my %POS2;
  foreach my $p (keys %BEST_READS_PER_POS2) {
     my $p2 = int($p/$window_for_merging)+1;
     push @{$POS2{$p2}},@{$BEST_READS_PER_POS2{$p}};
  }

  my %READS_USED;
  my %READS_NOT_MERGED;
  my $mergedn = 0;

  print "$current/$ctgs_to_process: $ctgid ($ctglen): there are ",scalar keys %POS1," and ",scalar keys %POS2," positions to merge\n";
  for my $reduce (['posF',\%POS1,'2'], ['posR',\%POS2,'1']) {
       my ($dir,$POS,$d) = @$reduce;
       mkdir $dir;
       chdir $dir;
       my @pos;
       if ($d == 2) { 
          @pos = sort {$a <=> $b} keys %{$POS};
       }
       else { 
          @pos = sort {$b <=> $a} keys %{$POS};
       }
       my $totpos = scalar @pos;
       my $pn = 0;
       print "$current/$ctgs_to_process: $ctgid ($ctglen): $dir $totpos\n";
       foreach my $p (@pos) {
          $pn++;
          my (%THISREADS,%THISQUAL,%THISCON,%READS_CONSIDERED);
          foreach my $fq (@{$POS->{$p}}) {
              my ($h1,$seq1,$hq1,$qual1,$h2,$seq2,$hq2,$qual2) = split(/\n/,$fq);
              $qual1 = convert_qual($qual1,0);
              $qual2 = convert_qual($qual2,0);
              my ($id1) = $h1=~/\@(\S+)/;
              $id1=~s/\/1/\.r/; #!!! ASSUME all libraries are stranded rf, swappimg names to have f/r as +/- strand
              $id1=~s/\/2/\.f/;
              my ($id2) = $h2=~/\@(\S+)/;
              $id2=~s/\/1/\.r/;
              $id2=~s/\/2/\.f/;
              my ($id) = $id1=~/(\S+)\./;
              my ($id_check) = $id2=~/(\S+)\./;
              unless (defined $id && defined $id_check && $id eq $id_check) {
                 die "not paired retrieval: $id $id_check";
              }
              $THISREADS{$id}{'1'}=">$id1\n$seq1\n";
              $THISREADS{$id}{'2'}=">$id2\n$seq2\n";
              $THISQUAL{$id}{'1'}=">$id1\n$qual1\n";
              $THISQUAL{$id}{'2'}=">$id2\n$qual2\n";
              my ($lib) = $h1=~/\@(L\d+)/;
              $THISCON{$id}=join(" ",$id1,$id2,$SIZES{$lib}{'min'},$SIZES{$lib}{'max'})."\n";
          }
          my $chunk = 1; #just carry-over from old code
          open OUT1,">$p\_c$chunk.1.fa";
          open QUAL1,">$p\_c$chunk.1.fa.qual";
          open OUT2,">$p\_c$chunk.2.fa";
          open QUAL2,">$p\_c$chunk.2.fa.qual";
          foreach my $id (keys %THISREADS) {
            print OUT1 $THISREADS{$id}{'1'};
            print OUT2 $THISREADS{$id}{'2'};
            print QUAL1 $THISQUAL{$id}{'1'};
            print QUAL2 $THISQUAL{$id}{'2'};
          }
          close OUT1; close OUT2; close QUAL1; close QUAL2;
          if (-s "$p\_c$chunk.1.fa" > 0) {
              print "$current/$ctgs_to_process: $ctgid ($ctglen): ",($p-1)*$window_for_merging,"-",$p*$window_for_merging," $dir $pn/$totpos,",scalar keys %THISREADS," reads\n";
              my ($r,$q,$con,$newpairn,@reads_used) = assemble_pairs("$p\_c$chunk",\%THISREADS,\%THISQUAL,$d,$unqid,$pairn);
              $pairn = $newpairn;
              if (@reads_used) {
                print OUT $r,$q,$con;
                print OUT join(" ",@reads_used),"\n";
                $mergedn++;
                foreach my $rid (@reads_used) {
                  $READS_USED{$rid}++;
                }
              }
              else {
                foreach my $id (keys %THISREADS) {
                  $READS_NOT_MERGED{$id}=$THISREADS{$id}{'1'}.$THISREADS{$id}{'2'}.$THISQUAL{$id}{'1'}.$THISQUAL{$id}{'2'}.$THISCON{$id};
                }
              }
          }
          `rm *`; #cleanup after each pos
       }
       chdir ".."; 
  }

  my $nonmergedn = 0;
  foreach my $id (keys %READS_NOT_MERGED) {
     print OUT $READS_NOT_MERGED{$id},$id,"\n" unless defined $READS_USED{$id};
     $nonmergedn++;
  }
  print OUT "\n";
  my $reads_after_merging = $nonmergedn + $mergedn;
  print STAT join("\t",$ctgid,$totreadsn,$nonmergedn,$mergedn,$reads_after_merging,sprintf("%.2f",$reads_after_merging/$totreadsn*100)."%"),"\n";
  `rm -r *`;
}


close BLOCK;
close OUT;
close STAT;

chdir $workdir;
`rm -R $tmpfsdir`;

print "\n\nALL DONE\n\n";


##########
sub convert_qual {
  my ($str,$rc) = @_;
  $str=~s/\n//;
  my @arr = split(//,$str);
  for my $i (0..$#arr) {
     $arr[$i]=ord($arr[$i])-33;
  }
  @arr = reverse(@arr) if $rc > 0;
#  my $sum = 0;
#  foreach my $q (@arr) { #not used anymore
#    $sum+=$q;
#  }
#  return (join(" ",@arr), $sum);
  return join(" ",@arr);
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


##########
sub assemble_pairs {
    my ($p,$READS,$QUAL,$d,$unqid,$pairn) = @_;
    `cap3 $p.$d.fa $cap3params >$p.$d.cap3`;
    open ACE,"$p.$d.fa.cap.ace" || die "no $p.$d.fa.cap.ace\n";
    local $/ = undef;
    my @ace = split(/\n\n+/,<ACE>);
    close ACE;
    local $/ ="\n";
    my $ctgname = "";
    my %CTGINFO;
    my @ctg_reads;
    foreach my $block (@ace) {
      if ($block=~/^CO (\S+) (\d+) (\d+)/) {
        $ctgname=$1;
        push @ctg_reads,[$ctgname,$3];
        $CTGINFO{$ctgname}{'CO'}=$block;
      }
      elsif ($block=~/^BQ/) {
        $CTGINFO{$ctgname}{'BQ'}=$block;
      }
      elsif ($block=~/^AF/) {
        $CTGINFO{$ctgname}{'AF'}=$block;
      }
    }
#    my $pairn = 0;
    my %REPORTED;
    my ($s,$q,$con)=("","","");
    return($s,$q,$con,$pairn,keys %REPORTED) unless @ctg_reads;
    my @highest_ctg = sort {$b->[1] <=> $a->[1]} @ctg_reads;
#   foreach my $ctg (keys %CTGINFO) {
    foreach my $ctg ($highest_ctg[0]->[0]) { #only highest coverage ctg is considered
      my $j = 0;
      open POUT1,">$p.$ctg.$j.1"; open PQUAL1,">$p.$ctg.$j.1.qual";
      open POUT2,">$p.$ctg.$j.2"; open PQUAL2,">$p.$ctg.$j.2.qual";
      while ($CTGINFO{$ctg}{'AF'}=~/AF (\S+)\.(\S) (\S) (\d+)/g) {
        my ($id,$direction,$orientation,$start) = ($1,$2,$3,$4);
        if ($orientation eq 'U') {
          print POUT1 $READS->{$id}{'1'}; print PQUAL1 $QUAL->{$id}{'1'};
          print POUT2 $READS->{$id}{'2'}; print PQUAL2 $QUAL->{$id}{'2'};
        }
      }
      close POUT1; close PQUAL1;
      close POUT2; close PQUAL2;
      while (1) {
        `cap3 $p.$ctg.$j.1 $cap3params -r 0  >$p.$ctg.$j.1.cap3`;
        `cap3 $p.$ctg.$j.2 $cap3params -r 0  >$p.$ctg.$j.2.cap3`;
         my ($ctg1,$qual1,$RD1,$remove1) = read_ace_Contig1("$p.$ctg.$j.1.cap.ace","r","U");
         my ($ctg2,$qual2,$RD2,$remove2) = read_ace_Contig1("$p.$ctg.$j.2.cap.ace","f","U");
         my %COMMONREADS;
         my $missed=0;
         my %FLANKS;
         foreach my $id (keys %{$RD1}) {
            if (defined $RD2->{$id}) { 
               $COMMONREADS{$id}++; 
               $FLANKS{$id}{'r'} = $RD1->{$id}; 
            } 
            else  { $missed=1; }
         }
         foreach my $id (keys %{$RD2}) {
            if (defined $RD1->{$id}) { 
               $COMMONREADS{$id}++;
               $FLANKS{$id}{'f'} = $RD2->{$id}; 
            } 
            else { $missed=1; }
         }
         if (scalar keys %COMMONREADS == 0) {
           last;
         }
         if ($missed > 0 || $remove1 > 0 || $remove2 > 0) {
            $j++;
            open POUT1,">$p.$ctg.$j.1"; open PQUAL1,">$p.$ctg.$j.1.qual";
            open POUT2,">$p.$ctg.$j.2"; open PQUAL2,">$p.$ctg.$j.2.qual";
            foreach my $id (keys %COMMONREADS) {
              print POUT1 $READS->{$id}{'1'}; print PQUAL1 $QUAL->{$id}{'1'};
              print POUT2 $READS->{$id}{'2'}; print PQUAL2 $QUAL->{$id}{'2'};
            }
            close POUT1; close PQUAL1;
            close POUT2; close PQUAL1;
         }
         else {
            $pairn++;
            my @minsize;
            my @maxsize;
            my $id="m-$unqid\_$pairn";
            $s.=">$id.r\n$ctg1\n>$id.f\n$ctg2\n";
            $q.=">$id.r\n$qual1\n>$id.f\n$qual2\n";
            foreach my $id (keys %COMMONREADS) {
              my ($lib) = $id=~/^(L\d+)/;
              push @minsize, $SIZES{$lib}{'min'} + $FLANKS{$id}{'r'} + $FLANKS{$id}{'f'} - 2;
              push @maxsize, $SIZES{$lib}{'max'} + $FLANKS{$id}{'r'} + $FLANKS{$id}{'f'} - 2;
              $REPORTED{$id}++;
            }
            @minsize = sort {$a <=> $b} @minsize;
            @maxsize = sort {$b <=> $a} @maxsize;
            $con.="$id.r $id.f $minsize[0] $maxsize[0]\n";
            last;
         }
      }
    }
    return($s,$q,$con,$pairn,keys %REPORTED);
}



##########
sub read_ace_Contig1 {
    my ($ace,$direction,$orinetation) = @_;
    open(ACE,$ace) || die "Can't read $ace\n";
    local $/ = undef;
    my @ace = split(/\n\n+/,<ACE>);
    close ACE;
    local $/ ="\n";
    my $thisctgname = "";
    my %THISCTG;
    foreach my $block (@ace) {
      if ($block=~/^CO (\S+)/) {
         $thisctgname=$1;
         $THISCTG{$thisctgname}{'CO'}=$block;
      }
      elsif ($block=~/^BQ/) {
        $THISCTG{$thisctgname}{'BQ'}=$block;
      }
      elsif ($block=~/^AF/) {
        $THISCTG{$thisctgname}{'AF'}=$block;
      }
    }
    my %RD;
    my $remove=0;
    my $ctgseq="";
    my $ctgqual="";
    if (defined $THISCTG{'Contig1'}) {
      $ctgseq=$THISCTG{'Contig1'}{'CO'};
      $ctgseq=~s/.*\n//; $ctgseq=~s/\n//g; $ctgseq=~s/\*//g;
      $ctgqual=$THISCTG{'Contig1'}{'BQ'};
      $ctgqual=~s/.*\n//; $ctgqual=~s/\s+/ /g;
      while ($THISCTG{'Contig1'}{'AF'}=~/AF (\S+)\.(\S) (\S) (\d+)/g) {
         my ($id,$direction,$orientation,$start) = ($1,$2,$3,$4);
         if ($direction eq $direction && $orientation eq $orientation) {
           $RD{$id}=$start;
         }
         else {
           $remove=1;
         }
       }
     }
     return ($ctgseq,$ctgqual,\%RD,$remove);
}
######


###########################################


##########
sub count_qual {
  my ($str) = @_;
  $str=~s/\n//;
  my @arr = split(//,$str);
  for my $i (0..$#arr) {
     $arr[$i]=ord($arr[$i])-33;
  }
  my $sum = 0;
  foreach my $q (@arr) { #not used anymore
    $sum+=$q;
  }
  return $sum;
}
##########

