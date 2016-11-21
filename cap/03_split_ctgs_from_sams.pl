#!/usr/bin/perl -w

use strict;
use Parallel::ForkManager;
use Config::General qw(ParseConfig);
my %CONF = ParseConfig("CONFIG");

my $K = $CONF{'HITS_TO_REPORT'};

my $ctgparts = $CONF{'CTGPARTS'};

my $pm = Parallel::ForkManager->new(7);

#mkdir "stat.03";

while (my $libdir = <LIBS_LS/*>) {
   $pm->start and next;
     my ($lib) = $libdir=~/LIBS_LS\/(\S+)/;
     my %CTG;
     my $n = 0;
     while (my $f = <$libdir/*.sams>) {
       $n++;
       warn "$n. $f\n";
       open F, $f;
       local $/ = "\n\n";
       while (my $seq = <F>) {
         $seq=~s/^(\S+)\n//;
         my $p = $1;
         $seq=~s/\n\n/\n/;
         push @{$CTG{$p}}, $seq;
       }
       close F;
       if ($n=~/[2468]00$/) {
         warn "writing down\n";
         foreach my $p (keys %CTG) {
          if (scalar @{$CTG{$p}} > 0) {
            open OUT,">>$libdir/$lib.$p.ctg";
            print OUT @{$CTG{$p}};
            close OUT;
            @{$CTG{$p}}=();
          }
         }
       }
#       `rm $f`;
     }
     warn "writing down final\n";
     foreach my $p (keys %CTG) {
       if (scalar @{$CTG{$p}} > 0) {
           open OUT,">>$libdir/$lib.$p.ctg";
           print OUT @{$CTG{$p}};
           close OUT;
       }
     }
   $pm->finish;
}

$pm->wait_all_children();


