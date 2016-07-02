#!/usr/bin/perl
use strict;
use warnings;

#################################################################################
# Written by Nick Kent, Jan 2011
# Fiddled with: Nick Kent, Dec 12th 2011
# Editied by Dan Pass 2015
#
# This script takes the .txt output files from SAMparser2.pl and calculates a
# frequency distribution for paired read insert size dyad position values within
# user-specified bins (defined in the variable $bin_width). It outputs the results
# as an .sgr file for rendering in IGB or converting into .wig for danpos.
#
################################################################################

use Getopt::Std;

my %options=();
getopts('i:o:b:p:a:h', \%options);

my $usage = "## USAGE: sgr_builder.pl -i infile.txt -o outfile.sgr -p chromosome_sizes.txt
[REQUIRED]
  -i infile.txt
  -o outfile.sgr
  -p {chromosome size file. Format: NAME[tab]size}
[OPTIONAL]
  -b {bin width (default: 10)}
  -a INT Do a bin average on the data (Default OFF, recommended 3).
  -h print this helpful help page

NOTE: Parameter file format:\nChr1  10002020\nChr2  1241414\nChr3 1308571\nABCDEFG  123456\n
NOTE 2: Make sure the chromosome names match in your reference otherwise everything fails!";

if (!%options || exists $options{h}){
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

print "~~~~~~~~~~~~~~~~~~~~~~~\n";

# Bin averaging
my $binavg = 1;
if(exists $options{a}){
  $binavg = $options{a};
}else{
    print "Not averaging over bins, remmember to pass e.g. -a 3 to average (or other value)\n";
}

my @chrsize;
if(exists $options{p}){
  print "Reading chromosome sizes\n";
  my $parafile = $options{p};
  open(PARA, '<', $parafile) or die "Unable to access Chromosome size parameter file: $parafile $!";
  chomp(@chrsize = <PARA>);
}else{
  print "Using Arabidopsis thaliana chromosomes as you didn't say otherwise! If it sharts throwing \"Uninitialised Value\" errors, check that chomosome names are the same as in the reference file.\n";
  my $parafile = "/home/sbi6dap/Projects/DansProcessingPipeline/scripts/core_scripts/Atha_chr_sizes.txt";
  print "Loading chromosome sizes from $parafile\n";
  open(PARA, '<', $parafile) or die "Unable to access Chromosome size parameter file: $parafile $!";
  chomp(@chrsize = <PARA>);
}

print "~~~~~~~~~~~~~~~~~~~~~~~\n";

################
# MAIN PROGRAM #
################

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
  # Modify value into bin (10bp by default)
  my $bin = (int($line[3]/$bin_width))*$bin_width;
  # Count into bins
  $binhash{$line[0]}{$bin} = $binhash{$line[0]}{$bin} + 1;
}
close(IN);

# Calculate the 3 bin moving average
if (exists $options{a}){
  open(OUT, '>', $outfile) or die "Unable to open $outfile: $!";
  # Export bin counts per chromosome;
  for my $chr (sort keys %binhash){
    print "Processing $chr\n";
    for my $binout (sort {$a <=> $b} keys %{$binhash{$chr}}){
        my $blast = $binout - 10;
        my $bnext = $binout + 10;
        my $binaverage;
        if (!exists($binhash{$chr}{$blast}) || !exists($binhash{$chr}{$bnext})){
            #print "cant find" . $binhash{$chr}{$binout-1} . "\n";
            next;
        }else{
            $binaverage = sprintf("%.3f",($binhash{$chr}{$blast} + $binhash{$chr}{$binout} + $binhash{$chr}{$bnext}) / $binavg);
        }
        print OUT "$chr\t$binout\t$binaverage\n";
    }
  }
}else{
  open(OUT, '>', $outfile) or die "Unable to open $outfile: $!";
  # Export bin counts per chromosome;
  for my $chr (sort keys %binhash){
    print "Processing $chr\n";
    for my $binout (sort {$a <=> $b} keys %{$binhash{$chr}}){
      print OUT "$chr\t$binout\t$binhash{$chr}{$binout}\n";
    }
  }
}
