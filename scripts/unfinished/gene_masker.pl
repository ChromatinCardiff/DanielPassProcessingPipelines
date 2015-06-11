#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Std;
my %options=();
getopts('i:o:m:x:f:', \%options);

my $usage = "gene_masker.pl -i input.txt [form: cols-positions, rows-genes] -o outfile.txt \
[optional] -g annotations.genepred [genepred annotations] ~~OR~~ -m genemask.txt, -f AT2G07732 [Force gene name, ignore first colum (useful if input is one gene from multiple samples)]\n";


open(IN, '<', $options{i}) or die "$usage:$!";
open(OUTFILE, '>', $options{o}) or die "$usage:$!";
my %maskhash;

if (exists $options{m}){
  open(MASK, '<', $options{m}) or die "$usage:$!";
  while(<MASK>){
    chomp;
    my @line = split('\t', $_);
    $maskhash{$line[0]} = $line[1];
  }
}

my $window = 5000;
my $offset = 1500;
my $binsize = 10;

if (exists $options{g}){
  open(SIZES, '<', $options{g}) or die "$usage:$!";
  while (<SIZES>){
    chomp;
    my @line = split('\t', $_);
    my @starts = split(',', $line[8]);
    my @ends = split(',', $line[9]);

    my $winc = 0;
    my $TSS = $starts[0];
    my $mask;
    #print "$line[0]\n";

    while ($winc <= $window){
      $winc++;
      my $einc = 0;
      my $hit = 0;

      while ($einc < scalar @starts){
        if (($winc + $starts[0] >= $starts[$einc]) && ($winc + $starts[0] <= $ends[$einc])){
          $hit = 1;
        }
        $einc++;
      }
      $mask .= $hit;
    }
    $maskhash{$line[0]} = $mask;
    print "$line[0]\t$mask\n";
  }
}

my $header = <IN>;
select(OUTFILE);
print $header;

while (<IN>){
  chomp;
  my @gene = split('\s', $_);
  my $genename;
  
  if (exists $options{f}){
    $genename = $options{f};
  }else{
    $genename = $gene[0]
  }

  if (exists $maskhash{$genename}){
    print "$gene[0]";
    print "\tNA" x ($offset /$binsize);

    my @pos = (split '', $maskhash{$genename});
    my $pinc = 0;
    while ($pinc <= $window){
      if ($pos[$pinc] == 1){
        #print "1";
        print "\t$gene[($pinc + $offset)/10]";
      }else{
        print "\tNA";
      }
      $pinc += $binsize;
    }
    print "\n";
  }else{
    print STDOUT "Gene $genename doesnt exist in exon mask file\n";
  }
}
