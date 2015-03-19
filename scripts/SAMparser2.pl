#!/usr/bin/perl
# Written by: Steff Adams, May 2010
# 
# Altered to read .sam files (change from original ELAND input) 8th June 2010: Nick Kent
# Altered to output only the .txt file 13th Sept 2010: Nick Kent
# Last fiddled with: Nick Kent, Aug 2011
# Added command flags, removed config file requirements, instructions on error & some defaults: Dan Pass, Mar 2015
#
# Usage: perl SAMparser.pl
#
# This script takes the chrn_info.txt files from chrgrep.sh and calculates the 
# centre positions of the paired reads for each chromosome within user-defined 
# size classes (+/- a user-defined window).
#
# The script takes configuration input from two sub files:
#
# 1. config.txt: modify this file to set in and out directory paths and the pwind
# variable. pwind = "particle window" - this is best left at 0.2 which means that
# the script will define particles +/- 20% of a specified size class. e.g 150bp 
# particles will be defined as having an ISIZE/end-to-end distance of between 120 
# and 180bp. You will notice that the pwind variable widens the size selection as
# the size class increases. This is important to encompass variation in linker length
# when dealing with poly-nucleosomes. The choice of 0.2 also prevents Part50, Part 100,
# Part150 and Part300 classes from overlapping. This is probably not important
# but prevents double counting within these functionally important size classes 
#
# 2. particle_size_range.txt: sets the particle sizes; simply input numbers delimited
# with commas.
#
# You will also need to set the $SAM_ID_flag variable as in SAMhistogram.pl
#
################################################################################
use strict; 
use warnings;
use Cwd;

##############################
# Get options

use Getopt::Std;
my %options=();
getopts('i:o:p:w:f:E', \%options);

my $usage = "## USAGE: SAMparser.pl [REQUIRED] -i infile.sam -o outdirectory -f {Flowcell ID} [OPTIONAL] -w {% window (default 0.2)} -p {particle sizes(comma separated, default 0,100,150,300,450)}";

if (!%options){
	print "~~\n$usage\n~~\n" and die;
}
################################################################################
# SET THE VARIABLES BELOW AS REQUIRED
################################################################################

my $infile = $options{i}; 
my $outdir = $options{o}; 
my $SAM_ID_flag = $options{f};

my @psizes;
my $pwind;

if(exists $options{w}){
	$pwind = $options{w}; 
}else{
	$pwind = 0.2;
}

if(exists $options{p}){
	@psizes = split(',', $options{p}); 
}else{
	@psizes = qw/0 100 150 300 450/;
}

# For SAM. Change to 142 & 284 respectively for ELAND
my $read_one_qual = 255;
my $read_two_qual = 0;


print "Infile: $infile\nOutdirectory: $outdir\nParticle size window (+/-): $pwind\nFlowcell ID: $SAM_ID_flag\n";
foreach (@psizes){
	print "Bins at: $_ bp\n";
}


################################################################################
# main program
###############################################################################

my $phigh;
my $plow;
my @handles;

# get the relevant info from the input file name
(my $prefix = $infile) =~ s/^\.sam//;

# print out some useful info
print ("\nFiltering '".$infile."' to particle size:\n");

# process input file for each partic size
# nick alters to make $pwind a percentage of the particle size
open(IN, "$infile") || die "Unable to open $infile: $!";
foreach my $psize (@psizes){
        open $handles[$psize], '>', "$outdir/$prefix.$psize" || die "Unable to open outfile: $!";
}

while(<IN>){
    # we want to ignore column headers
    if($_ =~ /^$SAM_ID_flag/){
        chomp($_);
        # split line by delimiter and store elements in an array
        my @line = split('\t',$_);
    
        foreach my $psize (@psizes){
        
            #print $psize."bp\n";
            $phigh = ($psize + ($pwind * $psize));
            $plow = ($psize - ($pwind * $psize));
            my $posfile = $prefix."_".$psize.".txt";
        
            #print "$line[4] $line[8] plow=$plow phigh=$phigh\n";
            if(($line[4] >= $read_one_qual)
               && (($line[8] > $plow && $line[8] < $phigh)
                   || ($line[8] < (-$plow) && $line[8] > (-$phigh)))){
                    
    			  # removed write out to filtered file (.fil)
               
                    # write out the dyad position data to the position file
                if($line[8] > 0){
                  my $pos = ($line[3] + ($line[8] * 0.5));
                  select $handles[$psize];
                  print "$line[2]\t$line[3]\t$line[8]\t$pos\n";
                }
            }
        }
    }    
}
#close(OUTPOS);
close(IN);
################################################################################


