#!/usr/bin/perl
# Written: Nick Kent,  Aug 2010
# Last updated: Nick Kent, 12 Sep 2014
# USAGE:- perl peakmarker.pl
# This script takes an .sgr file as an input, and calls peak centre/summit 
# bins above a single, but scalable, threshold.
# It then lists these bin positions with a y-axis value proportional to the scaled 
# read#frequency.
# The scaling value can be chosen to reflect differences in read depth between two
# experiments - either based on total depth or a SiteWriter-derived local depth.
# For simple peak counts, just leave at 1.00
#
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. 
#
# IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO: PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; DEGREE FLUNKING; LAUNCH OF STRATEGIC NUCLEAR WEAPONS; 
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE, 
# IGNORANCE, RANK INSANITY OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
################################################################################

use strict;
use warnings;
use Math::Round;

################################################################################
# SET THE VARIABLES BELOW AS REQUIRED
# $indir_path   - The directory containing the .sgr files to be processed
# $outdir_path  - The directory to store the .sgr peak output files
# $thresh  - The aligned read number noise threshold value
# $scale_factor - A proportion based on differences in read depth
################################################################################

my $indir_path ="/home/sbink/Nick_K/PeakMarker/sgr_in";
my $outdir_path ="/home/sbink/Nick_K/PeakMarker/peaks_out";
my $thresh = 30;
my $scale_factor = 1.0;


################################################################################
# MAIN PROGRAM
################################################################################

# define some variables
my (@files, $infile, $outfile, $n, @line);

# store input file names in an array
opendir(DIR, $indir_path) || die "Unable to access file at: $indir_path $!\n";
@files = readdir(DIR);

# process each input file within the indir_path in turn
foreach $infile (@files){
    
    # ignore hidden files and only get those ending .sgr
    if (($infile !~ /^\.+/) && ($infile =~ /.*\.sgr/)){
        
        # define outfile name from infile name
        $outfile = substr($infile,0,-4)."_peak_t".$thresh;
        $outfile .= '.sgr';
        
       
        # print out some useful info
        print ("\nProcessing '".$infile."'\n");
        
        open(IN, "$indir_path/$infile")
            || die "Unable to open $infile: $!";
        
        # define three new arrays to store required values from infile
        my @chr;
	my @bins;
        my @freq;
        
        # loop through infile to get values
        while(<IN>){
            
            # split line by delimiter and store elements in an array
            @line = split('\t',$_);
            # store the columns we want in two new arrays
            push(@chr,$line[0]);
	    push(@bins,$line[1]);
            push(@freq,$line[2]);
        }
        
        # close in file handle
        close(IN);
        
        # store size of array
        my $size = @bins;
        
        
        # try and open output file
        open(OUT,"> $outdir_path/$outfile")
             || die "Unable to open $outfile: $!";
        
        
        # need a variable to store line count
        my $count = 0;
        
        # this calls the peaks - giving an x-axis bin value ONLY for the peak centre and
	# a y-axis value as the peak hight scaled to some value proportionate to relative
	# read depth a relevant pair-wise comparison


	while ($count < $size){

	if (($freq[$count]>=$freq[$count-1]) && 
		($freq[$count]>$freq[$count+1]) && 
		($freq[$count]*$scale_factor>=$thresh)){
		
	  print(OUT 	$chr[$count]."\t".
			$bins[$count]."\t".
			round($freq[$count]*$scale_factor)."\n");
			
	$count++;
        

	}else{
	$count++;
        
        
            
        }}
        # close out file handle
        close(OUT);
    }
}