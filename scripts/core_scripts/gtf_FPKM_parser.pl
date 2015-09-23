#!/usr/bin/perl
use strict;
use warnings;

#################################################
# daniel.antony.pass@googlemail.com | July 2015 #
#################################################

##############################
# Get options
##############################
use Getopt::Std;
my %options=();
getopts('i:o:b:a:', \%options);
my $sep = "~~~\n";

my $usage = "$sep Take a GTF and sperate into multiple smaller gtf files based on FPKM values
USAGE: gtf_FPKM_parser.pl -i infile.gtf
[REQUIRED]
 -i [infile in gtf format]
[OPTIONAL]
 -o [outdir | default: .]
 -b [Percentage boundries of bins to split files into | default: 20,40,60,80 (therefore generating 5 files, <20%, 20%<x<40%, 40%<x<60%, 60%<x<80%, >80% )]
 -a [Absolute boundries of bins to split files into in FPKM | default: 1,10,100,1000 (therefore generating 5 files, <1%, 1%<x<10%, 10%<x<100%, 100%<x<1000%, >1000% ]\n";

if ((!%options) or ($options{h})){
	print "$usage$sep" and die;
}

# Infile name
my $infile = $options{i};
(my $prefix = $infile) =~ s/\.gtf//;

# Assign defaults
my $outdir = ".";
my @thresholds = qw/20 40 60 80 100/;

# Out directory
if(exists $options{o}){
	$outdir = $options{o};
	mkdir $outdir;
}

# thresholds
if(exists $options{b}){
	@thresholds = split(',', $options{b});
  push(@thresholds, "max");
  print "$sep" . "Splitting $infile into " . scalar @thresholds . " files based on percentage boundaries:\n";
}elsif(exists $options{a}){
  @thresholds = split(',', $options{a});
  push(@thresholds, "max");
  print "$sep" . "Splitting $infile into " . scalar @thresholds . " files based on nunmeric FPKM boundaries:\n";
}

foreach (@thresholds){
  print "$_\n";
}

### Main script

open(IN, "<", "$infile") || die "Unable to open $infile: $!";

# get total and boundary values
my ($linecount, @FPKMs);

while (my $line = <IN>){
  $linecount++;
  $line =~ /FPKM "([0-9]+.[0-9]+)";/;
  push(@FPKMs, $1);
}

my @sortedFPKMs = sort { $a <=> $b } @FPKMs;

my %FHs;
my %THs;

# Open outfiles ans determine threshold values
foreach my $th (@thresholds){
  my $filename = "$outdir/$prefix" . "-FPKM-lt-$th.gtf";
  open (my $fh, ">", $filename) or die "Unable to create $filename | $!\n";
  $FHs{$th} = $fh;

  # threshold values if working from percentage
  if ($th eq 'max'){
    $THs{$th} = $sortedFPKMs[-1];
  }else{
    if(exists $options{b}){
      $THs{$th} = $sortedFPKMs[$linecount * ($th / 100) - 1];
    }else{
      $THs{$th} = $th;
    }
  }
}

print "$sep" . "Splitting files at boundaries:\n";
for my $val (sort {$a <=> $b} values %THs){
  print "$val\n";
}

seek IN, 0, 0;

while (my $line = <IN>){
  $line =~ /FPKM "([0-9]+.[0-9]+)";/;
  my $FPKM = $1;
  #print $FPKM;

  foreach my $th (@thresholds){
    #print "$FPKM vs $THs{$th}\n";
    if ($FPKM < $THs{$th}){
      select $FHs{$th};
      print "$line";
      last;
    }
  }
}
