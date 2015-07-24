#!/usr/bin/perl
use strict;
use warnings;
use Math::Round;
##############################
# Get options
##############################
use Getopt::Std;
my %options=();
getopts('i:b:o:h', \%options);
my $sep = "~~~\n";
my $usage = "## USAGE:  size_profiler.pl -i  size_at_positions.sgr
  [REQUIRED]
    -i infile.sgr (columns: 1-chromosome 2-position 3-read size)
  [OPTIONAL]
    -o outfile.sgr
    -b binsize for averaging (Default: 10)\n";


my %binhash;
my %counthash;
my $binsize = 10;
my $infile;

if (exists $options{h}){
  print $usage;
}

if (exists $options{i}){
  open(IN, "<", $options{i}) or die "Unable to open infile: $usage";
}else{
  print "No infile provided!\n$usage";
  die;
}

if (exists $options{o}){
  open(OUT, ">", $options{o}) or die "Unable to open outfile: $usage";
  select OUT;
}

# Bin size
if(exists $options{b}){
	$binsize = $options{b};
}

while(my $line = <IN>){
  # split line by delimiter and store elements in an array
  my @line = split('\t',$line);
  my $chr = $line[0];
  my $bin = (int($line[1]/10))*10;
  my $size = $line[2];

  if (exists $binhash{$chr}{$bin}){
    $binhash{$chr}{$bin} += $size;
  }else{
    $binhash{$chr}{$bin} = $size;
  }
  $counthash{$chr}{$bin}++;
}

for my $chr (sort keys %binhash){

  my $last;
  for my $bin (reverse sort {$a<=>$b} keys %{ $binhash{$chr} }){
    $last = $bin;
    #print "last is $last\n";
    last;
  }

  for (my $inc = 0; $inc <= $last; $inc += $binsize){
    #print "$inc\n";
    my $binval = 0;

    if (exists $binhash{$chr}{$inc}){
      $binval = round($binhash{$chr}{$inc} / $counthash{$chr}{$inc});
    }
    print  "$chr\t$inc\t$binval\n";
  }
}
