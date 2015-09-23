#!/usr/bin/perl
use strict;
use warnings;

##################################################
## Daniel Pass | github.com/passdan | August 2015 ##
##################################################

###############
# Get options #
###############

use Getopt::Std;
my %options=();
getopts('i:o:a:b:X', \%options);
my $sep = "~~~\n";
my $usage = "## USAGE:  plus1_slider.pl -i feature_positions.txt
  [REQUIRED]
    -i infile.txt (columns: samples, Rows: positions)
  [OPTIONAL]
    -o outfile.txt
    -a open window for slide detection (Default: 150 (Zero in 1.5kb in 10bp bins))
    -b close window for slide detection (Default: 170 (200 in 1.5kb in 10bp bins))
    -X output the value of the shift to the peak only, not the actual trace\n";

my ($a,$b,$root);

if (exists $options{i}){
  open(IN, "<", $options{i}) or die "Unable to open infile: $usage";
}else{
  print "No infile provided!\n$usage";
  die;
}

if (exists $options{o}){
  $root = $options{o};
}else{
  $root = $options{i} . "-1nucl_align";
}

if (exists $options{X}){
  open(SHIFTOUT, ">", "$root-shift.txt") or die "Unable to open outfile: $usage";
}else{
  open(OUT, ">", "$root-full.txt") or die "Unable to open outfile: $usage";
  open(AVGOUT, ">", "$root-avg.txt") or die "Unable to open outfile: $usage";
}

if (exists $options{a}){
  $a = $options{a};
}else{
  print "Default window opening used: 150\n";
  $a = 150;
}

if (exists $options{b}){
  $b = $options{b};
}else{
  print "Default window closing used: 170\n";
  $b = 170;
}

#############
# Main body #
#############

my $windowsize = $b - $a;
#select OUT;

my %posavg; #averaging hash
my $lc; #Line count

my $header = <IN>;
my @header = split(' ', $header);
my @h = @header[0 .. ($#header - $windowsize)];
print OUT join ("\t", @h) . "\n";

while(my $line = <IN>){
  $lc++;
  my @line = split(' ', $line);
  my $annot = shift @line;
  my @window = @line[$a...$b];

  my $i = 0;
  my $max =0;
  my $maxpos = 0;

  foreach my $bin (@window){
    if ($max < $bin){
      $max = $bin;
      $maxpos = $i;
    }
    $i++;
  }
  #print "Max was $max at position $maxpos\n";
  my @line_fix = @line[$maxpos .. ($#line -($windowsize - $maxpos))];

  my $inc=0;
  foreach my $p (@line_fix){
    $posavg{$inc} += $p;
    #print "$p";
    $inc++;
  }

  #splice @line, , $maxpos;

  print "$annot\t" . join ("\t", @line_fix) . "\n";
  #print "$annot\t" . "NA\t" x ($windowsize - $maxpos) . join ("\t", @line) . "\tNA" x $maxpos . "\n";

}
close OUT;

select AVGOUT;
print join ("\t", @h) . "\n";
print "$options{i}\t";

foreach my $val (sort keys %posavg){
  print $posavg{$val} / $lc . "\t";
}
