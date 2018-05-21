#!/usr/bin/perl
use strict;
use warnings;

##################################################
## Daniel Pass | github.com/passdan | August 2015 ##
##################################################

require Math::Spline;

###############
# Get options #
###############

use Getopt::Std;
my %options=();
getopts('i:o:a:b:s:X', \%options);
my $sep = "~~~\n";
my $usage = "## USAGE:  spline_slider.pl -i feature_positions.txt
  [REQUIRED]
    -i infile.txt (columns: samples, Rows: positions)
  [OPTIONAL]
    -o output-handle | The root for output files (output-handle-full.txt, output-handle-avg.txt, output-handle-shift.txt. Default: infile-1nucl)
    -a INT | Open bin position for window of slide detection (Default: 150 (TSS for 1.5kb in 10bp bins))
    -w INT | Window size for slide detection (Default: 20 (200bp in 10bp bins))
    -s shiftvalues.txt | Use a file generated with -X to decide on the shift, do not calulate again
    -X output the value of the shift to the peak only, not the actual trace\n";

my ($a,$b,$root,$windowsize);

# Input file
if (exists $options{i}){
  open(IN, "<", $options{i}) or die "Unable to open infile: $usage";
}else{
  die "$sep No infile provided!\n$sep$usage";
  #die;
}

# If using a precalculated shift file
my @shifts;
if (exists $options{s}){
  open(SHIFTIN, "<", $options{s}) or die "Unable to open reference shift file: $usage";
  @shifts = <SHIFTIN>;
}

# Output files
if (exists $options{o}){
  $root = $options{o};
}else{
  $root = $options{i} . "-1nucl_align";
}

# If you intend to only output a shift file for future funzies
if (exists $options{X}){
  open(SHIFTOUT, ">", "$root-shift.txt") or die "Unable to open outfile: $usage";
}else{
  open(OUT, ">", "$root-full.txt") or die "Unable to open outfile: $usage";
  open(AVGOUT, ">", "$root-avg.txt") or die "Unable to open outfile: $usage";
}

# For manual window changing
if (exists $options{a}){
  $a = $options{a};
}else{
  print "Default window opening position used: 150 bins\n";
  $a = 150;

}
# For manual window changing
if (exists $options{w}){
  $windowsize = $options{w};
}else{
  print "Default window size used: 20 bins\n";
  $windowsize = 20;
}

#############
# Main body #
#############

#select OUT;

my %posavg; #averaging hash
my $lc; #Line count
my (@header,@h); #header

if (!exists $options{X}){       # only print header if not doing the 'shift only'
  my $header = <IN>;
  @header = split(' ', $header);
  push(@h, shift(@header));

  foreach my $x (@header){
    $x = $x * 10;
  }

  my $last = (@header[($#header + 1 - $windowsize)]);
  for (my $i = $header[0] ; $i <= $last; $i++) {
    push (@h, $i);
  }
  print OUT join ("\t", @h) . "\n";
}

while(my $line = <IN>){
  $lc++;
  my @line = split(' ', $line);
  my $annot = shift @line;
  #print STDOUT @line;
  #print STDOUT @header;
  my @spline;
  my $spline = Math::Spline->new(\@header,\@line);
  my $count =0;
  my $position = ($header[0]);
  while ($position <= (($header[-1]))){
    #if($spline->evaluate($position) >0){
      push @spline, $spline->evaluate($position);
  #  }else{
  #    push @spline, "0";
  #  }

    #print "$position: @spline[$position]\n";
    $count++;
    $position++;
  }

  #print $spline->evaluate(-1551) . "\t";

  my @window = splice(@spline, ($a * 10),($windowsize * 10));    # extract window of investigation x spline expansion

  #print scalar @line . "\t";
  #print scalar @spline . "\t";
  #print scalar @window . "\t";
  #print $count . "\n";

  my $maxpos = 0;   ##The important number

  # use Max position as previously outputted by -X option
  if (exists $options{s}){
    my $preshifts = $shifts[$lc];
    $maxpos = (split /\t/, $preshifts)[1];

    ## Finding maximum position within defined range
  }else{
    my $i = 0;
    my $max = 0;

    foreach my $bin (@window){
      if ($max < $bin){
        $max = $bin;    # If this bin is bigger than the current max, then replace
        $maxpos = $i;   # If this bin is bigger than the current max, then assign $maxpos to this bin number
      }
      $i++;
    }
  }

  if (exists $options{X}){
    print SHIFTOUT "$annot\t$maxpos\n";
  }else{
    # Cut $maxpos bins from the left side, and the remainder of the $windowsize from the right
    my @line_fix = splice(@spline, $maxpos, ($#spline + 1 - $windowsize + $maxpos));

    my $inc=0;
    foreach my $p (@line_fix){
      $posavg{$inc} += $p;
      #print "$p";
      $inc++;
    }

    print OUT "$annot\t" . join ("\t", @line_fix) . "\n";
    #print "$annot\t" . "NA\t" x ($windowsize - $maxpos) . join ("\t", @line) . "\tNA" x $maxpos . "\n";
  }
}

if (!exists $options{X}){       # only print header if not doing the 'shift only'
  select AVGOUT;
  print join ("\t", @h) . "\n";
  print "$options{i}\t";

  foreach my $val (sort keys %posavg){
    print $posavg{$val} / $lc . "\t";
  }
}
