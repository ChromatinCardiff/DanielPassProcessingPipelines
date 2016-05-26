#!/usr/bin/perl
# Written by: Steff Adams, May 2010
#
# Nick Kent
# Altered to read .sam files (change from original ELAND input) 8th June 2010: Nick Kent
# Altered to output only the .txt file 13th Sept 2010: Nick Kent
#
# Dan Pass, Mar 2015
# Improved speed, added command flags, removed config file requirements, instructions on error, some defaults:
#
################################################################################
use strict;
use warnings;
use Cwd;

##############################
# Get options
##############################
use Getopt::Std;
my %options=();
getopts('i:o:p:w:f:E', \%options);

my $usage = "## USAGE: Example: SAMparser.pl -i infile.sam -o output_directory -f HISEQ
	[REQUIRED]
		-i infile.sam
		-o outdirectory
		-f {Flowcell ID}
	[OPTIONAL]
		-w {window size by percentage (No default, suggested 0.2)}
		-W {window size in absolute bp (default: 10)}
		-p {particle sizes(comma separated, default 100,150,300,450)}
		-E {will change from SAM quality (default) to ELAND}\n";

if (!%options){
	print "~~\n$usage\n~~\n" and die;
}

################################################################################
# Parameter set up
################################################################################

my $infile = $options{i};
my $SAM_ID_flag = $options{f};
my $outdir;
my @psizes;
my $pwind;

# get the relevant info from the input file name
(my $p1 = $infile) =~ s/^(.*\/)(\w+\.sam)$/$2/;
my $def_outdir = $1;
(my $prefix = $p1) =~ s/^(\w+)\.sam$/$1/;

# Assign default values
if(exists $options{o}){
	$outdir = $options{o};
}else{
	$outdir = "$def_outdir" . "particles";
}

# Make output directory if doesnt exist
if (-e $outdir){
	print "$outdir exists\n";
}else{
	my $cwd = getcwd;
	print "making directory $cwd/$outdir\n";
	mkdir "$cwd/$outdir";
}

if(exists $options{w}){
	$pwind = $options{w};	# Percentage size window #
}elsif(exists $options{W}){
	$pwind = $options{W}	# Fixed size window #
}else{
	$pwind = 10;
}

if(exists $options{p}){
	@psizes = split(',', $options{p});
}else{
	@psizes = qw/100 150 300 450/;
}

# Default for for SAM. If using ELAND, pass the -E flag
my $read_one_qual = 255;
my $read_two_qual = 0;

if(exists $options{E}){
	$read_one_qual = 142;
	$read_two_qual = 284;
}

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


# print out some useful info
print ("\nFiltering '".$infile."' to particle size:\n");

# process input file for each partic size
# nick alters to make $pwind a percentage of the particle size
open(IN, "$infile") || die "Unable to open $infile: $!";
foreach my $psize (@psizes){
	print "Opening: $outdir/$prefix" . "_$psize.txt\n";
        open ($handles[$psize], '>', "$outdir/$prefix" . "_$psize.txt") or die "Unable to open outfile: $!";
}

while(<IN>){
    # we want to ignore column headers
    if($_ =~ /^$SAM_ID_flag/){
        chomp($_);
        # split line by delimiter and store elements in an array
        my @line = split('\t',$_);

        foreach my $psize (@psizes){

            #print $psize."bp\n";
			if(exists $options{w}){
            	$phigh = ($psize + ($pwind * $psize));
            	$plow = ($psize - ($pwind * $psize));
			}else{
				$phigh = ($psize + $pwind);
				$plow = ($psize - $pwind);
			}
            my $posfile = $prefix."_".$psize.".txt";

            #print "$line[4] $line[8] plow=$plow phigh=$phigh\n";
            if(($line[4] >= $read_one_qual)
               && (($line[8] > $plow && $line[8] < $phigh)
                   || ($line[8] < (-$plow) && $line[8] > (-$phigh)))){

                # write out the dyad position data to the position file
                if($line[8] > 0){
                  my $pos = ($line[3] + ($line[8] * 0.5));
                  select $handles[$psize];

				#Chromosome | dyad start | size | apex position
                  print "$line[2]\t$line[3]\t$line[8]\t$pos\n";
                }
            }
        }
    }
}
close(IN);
