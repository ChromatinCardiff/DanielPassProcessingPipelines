#!/usr/bin/perl
# Written by Nick Kent, Jan 2011
# Fiddled with: Nick Kent, Dec 12th 2011
# USAGE:- perl Athal_histogramn.pl
#
# This script takes the .txt output files from nksamparserpro.pl and calculates a
# frequency distribution for paired read insert size dyad position values within
# user-specified bins (defined in the variable $bin_width). It outputs the results
# as an .sgr file for rendering in IGB.
#
# This script replaces nkhistogram.pl which was just soooooooo slow with big files.

# The implementation here is just a simple tally counter. The section starting at
# Line 57 has been modified (from the original Dicty version) to describe the
# lengths of each yeast chromosome.
#
# Clunky programming, BUT it works!
################################################################################
use strict;
use warnings;
#use Math::Round;
use Getopt::Std;

my %options=();
getopts('i:o:b:p:', \%options);

my $usage = "## USAGE: sgr_builder.pl [REQUIRED] -i infile.txt -o outfile.sgr -p {chromosome size file. Format: NAME[tab]size} -b {bin width (default: 10)}\n
              Parameter file format:\nChr1  10002020\nChr2  1241414\nChr3 1308571\nABCDEFG  123456\n";

if (!%options){
  print "~~\n$usage\n~~\n" and die;
}

my $infile = $options{i};
my $outfile = $options{o};

# Bin width default
my $bin_width = 10;

# Replace bin_width
if(exists $options{b}){
  my $bin_width = $options{b};
}

my @chrsize;
if(exists $options{p}){
  print "Reading chromosome sizes\n";
  my $parafile = $options{p};
  open(PARA, '<', $parafile) or die "Unable to access Chromosome size parameter file: $parafile $!";
  chomp(@chrsize = <PARA>);
}

################################################################################
# MAIN PROGRAM - you should not need to edit below this line
################################################################################

# define some variables
my ($n, @line);

# store input file names in an array
open(IN, '<', $infile) or die "Unable to access input file: $infile $!\n";

my %binhash;

# Populate hash with 0s upto max for each chromosome
foreach my $line (@chrsize){

  my ($name,$chrmax) = split('\t',$line);
  print "Populating Hash for $name at size $chrmax\n";

  for (my $i = 0; $i <= ($chrmax + $bin_width); $i+= $bin_width){
    $binhash{$name}{$i} = 0;
    #print "name: $name\tbin: $i\n";
  }
}


print "Calculating bin abundances\n";
while(<IN>){
  # split line by delimiter and store elements in an array
  @line = split('\t',$_);
  my $bin = (int($line[3]/$bin_width))*$bin_width;
  #print "Chr: $line[0]\tbin:$bin\tPos:$line[3]\n";
  $binhash{$line[0]}{$bin} = $binhash{$line[0]}{$bin} + 1;
}
close(IN);

open(OUT, '>', $outfile) or die "Unable to open $outfile: $!";
# Export bin counts per chromosome;
for my $chr (sort keys %binhash){
  print "Processing $chr\n";
  for my $binout (sort {$a <=> $b} keys %{$binhash{$chr}}){
    print OUT "$chr\t$binout\t$binhash{$chr}{$binout}\n";
  }
}
