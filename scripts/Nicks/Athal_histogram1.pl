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
use Math::Round;
################################################################################
# SET THE 3 VARIABLES BELOW AS REQUIRED
# $indir_path   - The directory containing the .txt files to be processed
# $outdir_path  - The directory to store the .sgr output files
# $bin_width    - The bin width value
################################################################################

my $indir_path ="/home/sbink/Nick_K/Histogram/histo1/in";
my $outdir_path ="/home/sbink/Nick_K/Histogram/histo1/out";
my $bin_width = 10;

################################################################################
# MAIN PROGRAM - you should not need to edit below this line
################################################################################

# define some variables
my (@files, $infile, $outfile, $n, @line);

# store input file names in an array
opendir(DIR, $indir_path) || die "Unable to access file at: $indir_path $!\n";
@files = readdir(DIR);

# process each input file within the indir_path in turn
foreach $infile (@files){
    
    # ignore hidden files and only get those ending .txt
    if (($infile !~ /^\.+/) && ($infile =~ /.*\.txt/)){
        
        # define outfile name from infile name
        $outfile = substr($infile,0,-4)."_".$bin_width;
        $outfile .= '.sgr';
        
        # get chromosome number from infile name
        $infile =~ /\w{3}(\d+).*/;
        $n = $1;
		
		# a cumbersome way to define bin maxima for TAIR_10 Chr 1-5 and ChrC and M
		my $top = 0;
		if ($n == 1){
			$top = 30427671;
			}
			elsif ($n==2){
			$top = 19698289;
			}
			elsif ($n==3){
			$top = 23459830;
			}
			elsif ($n==4){
			$top = 18585056;
			}
			elsif ($n==5){
			$top = 26975502;
			}
			elsif ($n == 6){
			$top = 154478; #ChrC
			}
			elsif ($n==7){
			$top = 366924; #ChrM
			}
			else {
			print "This is not Arabidopsis. There are too many chromosomes\n";
			exit;
			}
			
        # print out some useful info
        print ("\nProcessing '".$infile."'\n");
        
        open(IN, "$indir_path/$infile")
            || die "Unable to open $infile: $!";
        
        # define new array to store required dyad position values from infile
        my @dyad_pos;
        
        # loop through infile to get values
        while(<IN>){
            
            # split line by delimiter and store elements in an array
            @line = split('\t',$_);
            # store the column we want in a new array
            push(@dyad_pos,$line[3]);
        }
        
        # close in file handle
        close(IN);
        # Nick changes from here
		
		my $dyadarray_size= @dyad_pos;
		
		# Define the number of bins for the relevant chromosome
		my$bin_no = (int($top/$bin_width))+1;
		
		# Define the distribution frequency array
		my @dist_freq;
		my $i=0;
		
		# Fill the frequency distribution "bins" with zeroes
		for ($i=0; $i<$bin_no; $i++){
			push (@dist_freq, 0);
			}
			
		# Reset the incrementor and define the bin hit variable
		$i=0;
		my $bin_hit = 0;
		
		# The tally counter 
		while ($i < $dyadarray_size){
			$bin_hit = int($dyad_pos[$i]/$bin_width);
			$dist_freq[$bin_hit] ++;
			$i ++;
			}
		
		# Calculate the 3 bin moving average
		my @moving_ave;
		my $ma = 0;
		my $count = 1;
		push (@moving_ave,0);
		
		while ($count<$bin_no-1){
			$ma = (($dist_freq[$count--] + $dist_freq[$count] + $dist_freq[$count++])/3);
			push (@moving_ave,$ma);
			$count ++;
			}
			push (@moving_ave,0);
			
		
		
        # try and open output file
        open(OUT,"> $outdir_path/$outfile")
             || die "Unable to open $outfile: $!";
        
        
            # print required data in tab-delimited format to output file
            # NK modifies to output chrn, bin and ma only
			for ($i=0; $i<$bin_no; $i++){
			
            print(OUT "chr".$n."\t".($i*$bin_width)."\t".round($moving_ave[$i])."\n");
			}
            
        }
        # close out file handle
        close(OUT);
    
}