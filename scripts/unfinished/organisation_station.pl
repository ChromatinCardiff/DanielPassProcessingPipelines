#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Std;
my %options=();
getopts('i:o:m:x:b:', \%options);

my $usage = "organisation_station.pl -i input.txt [form: cols-positions, rows-genes] -o outfile.txt -m genemask.txt -x offset\n";

open(IN, '<', $options{i}) or die "$usage:$!";
open(DAVGOUT, '>', "$options{o}_avg.txt") or die "$usage:$!";
open(DLSTOUT, '>', "$options{o}_list.txt") or die "$usage:$!";
open(MASK, '<', $options{m}) or die "$usage:$!";

my ($offset, $binsize);

if(exists $options{x}){
  $offset = $options{x};
}else{
  $offset = 1500;
}

if(exists $options{b}){
  $binsize = $options{b};
}else{
  $binsize = 10;
}

my %maskhash;

print STDOUT "Importing genemask\n";
while(<MASK>){
  chomp;
  my @line = split('\t', $_);
  $maskhash{$line[0]} = $line[1];
}

my $head = <IN>;
#print $head;

#hashes for ordering
my %dhash_avg;
my %dhash_values;

#Incrementer for showing progress
my $sinc;
print STDOUT "Processing input\n";

while (<IN>){
  chomp;
  $sinc++;

  my @gene = split('\s', $_);
  if (exists $maskhash{$gene[0]}){
    #print "$gene[0]";

    if (($sinc / 1000) =~ /^\d+$/){
      print STDOUT "Processing line $sinc...\n";
    }

    my @maskpos = (split '', $maskhash{$gene[0]});

    #Position increment
    my $pinc = 0;
    my ($gsum, $glen);
    my @exonheights;

    #calculate gene mean in exonic regions (or median)
    foreach my $maskpos (@maskpos){
      $pinc++;
      #get last '1' (aka final exon position) for gene length
      if ($maskpos == 1){
        $glen = $pinc;
      }
      #print "$pinc $binsize " . $pinc/$binsize . "\n";

      if (($pinc / $binsize) =~ /^\d+$/){
        if ($maskpos == 1){
          my $genepos = (($pinc + $offset) / $binsize);
          push(@exonheights, $gene[$genepos]);      ## FOR MEDIAN CALCULATION
          #$gsum += $gene[$genepos];                ## FOR MEAN CALCULATION
        }
      }
    }
    # ERROR: Fails if gene length is less than 200ish DONT KNOW WHY
    if ($glen < 200){
      next;
    }
    my $gmean = $gsum / $glen;

    # Median calculation
    @sortexonheights = sort {$a <=> $b} @exonheights;
    my $gmedian = @sortexonheights[]
    #print "$gene[0]\tgmean: $gmean\tgsum: $gsum\tglen: $glen\n";

    #Calculate d at each position
    my $dgen = $gene[0];
    my $dsum;
    $pinc = 0;

    foreach my $maskpos (@maskpos){
      $pinc++;
      my $d;
      if (($pinc / $binsize) =~ /^\d+$/){
        if ($maskpos == 1){
          my $genepos = (($pinc + $offset) / $binsize);
          #print "pinc: $pinc \t offset: $offset\n";
          #print "$pinc:\t abs(($gene[int(($pinc + ($offset / $binsize)) / $binsize)] - $gmean) / $gmean)\n";
          $d = abs(($gene[$genepos] - $gmean));# / $gmean);      ### EDIT divided by mean removed
          $dsum += $d;
        }else{
          $d = "NA";
        }
        #print "$d\t";
        $dgen .= "\t$d";
      }
    }
    #print "\n";
    $dhash_avg{$gene[0]} = ($dsum / $glen);
    print DLSTOUT "$dgen\n";
    #print "Average d value for $gene[0] is " . $dsum / $glen . "\n";
  }
}

print STDOUT "Sorting dispersion per gene average\n";
foreach my $key (reverse sort { $dhash_avg{$a} <=> $dhash_avg{$b} } keys %dhash_avg){
  print DAVGOUT "$key\t$dhash_avg{$key}\n";
}
close(DAVGOUT);
close(DLSTOUT);
