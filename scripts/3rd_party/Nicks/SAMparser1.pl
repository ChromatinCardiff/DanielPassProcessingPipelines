#!/usr/bin/perl
# Written by: Steff Adams, May 2010
# 
# Altered to read .sam files (change from original ELAND input) 8th June 2010: Nick Kent
# Altered to output only the .txt file 13th Sept 2010: Nick Kent
# Last fiddled with: Nick Kent, Aug 2011
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

################################################################################
# SET THE VARIABLES BELOW AS REQUIRED

# $SAM_ID_flag 	- The Illumina Flowcell ID number, or a unique substring
################################################################################


my $SAM_ID_flag = "HWI-D00";


################################################################################
sub load_param($); # returns a parameter value from the config file
sub load_psizes(); # returns user-specified "particle sizes" classes
################################################################################

# define parameter config file
my $paramfile = "config.txt";
# define partical size range file
my $psizefile = "partic_size_range.txt";
# define current working directory
my $cwd = getcwd;

################################################################################
# function to return user-specified parameters from the config file
################################################################################

sub load_param($){

    my $param = shift;   
    my $value;
    unless (-f $paramfile){
        print localtime().": Config file '".$paramfile."' ";
        print "cannot be found at $cwd\n";
        exit 1;
    }
    open(IN, $paramfile) || die "Unable to open $paramfile: $!";
    while(<IN>){
        if(/^\s*([^=]+)\s*=+\s*([^=]+)\s*$/){            
            if($param eq $1){
                $value = $2;
                chomp($value);
                return $value;
            }
        }
    }
    close($paramfile);
    
}
################################################################################
# function to return user-specified values from particle size range file 
################################################################################

sub load_psizes(){

    my @values;
    unless (-f $psizefile){
        print localtime().": Particle sizes file '".$psizefile."' ";
        print "cannot be found at $cwd\n";
        exit 2;
    }
    open(IN, $psizefile) || die "Unable to open $psizefile: $!";
    while(<IN>){
        if(/^\d+.*/){
            @values = split(',',"$_");
        }
    }
    close($psizefile);
    return @values;
    
}
################################################################################
# main program
###############################################################################

# define vars
my ($phigh,$plow,$psize,$posfile,$pos,@line,@files,$infile);
my $indir_path = load_param("indir_path");
my $outdir_path = load_param("outdir_path");
my $read_one_qual = load_param("read_one_qual");
my $read_two_qual = load_param("read_two_qual");
my $pwind = load_param("pwind");
my @psizes = load_psizes(); 

# store input file names in an array
opendir(DIR, $indir_path) || die "Unable to access file at: $indir_path $!\n";
@files = readdir(DIR);

# process each input file in the indir_path in turn
foreach $infile (@files){
    
    # ignore hidden files
    if ($infile !~ /^\.+/){

        # get the relevant info from the input file name
        $infile =~ /^\w{3}(\d+)\_(\d+)\_(.*)\.{1}\w{3}$/;
        my $n = "chr".$1;
        my $x = $2;
        my $strain = $3;

        # print out some useful info
        print ("\nFiltering '".$infile."' to particle size:\n");
        
        # process input file for each partic size
        # nick alters to make $pwind a percentage of the particle size
        foreach $psize (@psizes){
            
            print $psize."bp\n";
            $phigh = ($psize + ($pwind * $psize));
            $plow = ($psize - ($pwind * $psize));
            $posfile = $n."_".$x."_Part".$psize."_".$3.".txt";
            
            # open input/output files for reading and writing
            open(IN, "$indir_path/$infile")
                || die "Unable to open $infile: $!";
            open(OUTPOS,"> $outdir_path/$posfile")
                || die "Unable to open $posfile: $!";
            
            # loop though input file line by line
            # nick altered here to account for SAM read quality
            while(<IN>){
                # we want to ignore column headers
                if($_ =~ /^$SAM_ID_flag/){
                    chomp($_);
                    # split line by delimiter and store elements in an array
                    @line = split('\t',$_);
                    
                    #print "$line[4] $line[8] plow=$plow phigh=$phigh\n";
                    if(($line[4] >= $read_one_qual)
                       && (($line[8] > $plow && $line[8] < $phigh)
                           || ($line[8] < (-$plow) && $line[8] > (-$phigh)))){
                        
			  # removed write out to filtered file (.fil)
                    
                   
                        # write out the dyad position data to the position file
                        if($line[8] > 0)
                              {
                              $pos = ($line[3] + ($line[8] * 0.5));
                              print (OUTPOS "$n\t$line[3]\t$line[8]\t$pos\n");
                        }
                       
                    }
                    
                }  
            }
            
            # close file handlers
            close(OUTPOS);
            close(IN);
            
        }
    }
}
################################################################################


