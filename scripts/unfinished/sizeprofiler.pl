#!/usr/bin/perl

# Daniel Pass daniel.antony.pass@googlemail.com 16/04/2015

use strict;
use warnings;

##############################
# Get options
##############################
use Getopt::Std;
my %options=();
getopts('i:o:m:', \%options);

my $usage = "## USAGE: sizeprofiler.pl [REQUIRED] -i infile.txt -o outfile.txt -m [minimum gap] /n";

if (!%options){
	print "~~\n$usage\n~~\n" and die;
}

my $infile = $options{i};
my $outfile = $options{o};
my $mingap = $options{m};

open(IN, "$infile") || die "Unable to open $infile: $!";
open(OUT, >, "$outfile") || die "Unable to open $outfile: $!";

my $previous = 0;
my $current = 0;

select(OUT);
while (my $line = <IN>){
  $current = (split '\t',$line)[1];
  my $gap = $current - $previous;
  if ($gap > $mingap){
    chomp($line);
    print "$line\t$gap\n";
  }

  $previous = (split '\t',$line)[1];
}
