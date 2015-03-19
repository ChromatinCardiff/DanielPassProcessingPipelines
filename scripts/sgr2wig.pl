#!/usr/bin/perl 
use strict;
use warnings;

my $usage = "## USAGE:  sgr2wig.pl inputfile.sgr outputfile.wig\n";

my $infile = $ARGV[0];
my $outfile = $ARGV[1];

open(SGRIN, "$infile") or die $usage, $!;
open(WIGOUT, ">", "$outfile") or die $usage, $!; 

my $chr="chr"; 

select(WIGOUT);
while(<SGRIN>){
        my @line=split/\s+/, $_; 

        if ($line[0] ne $chr){ 
		# New chromosome line
		print "fixedStep chrom=$line[0] start=0 step=10 span=10\n";
		$chr = $line[0];
		print STDOUT "## Processing: $line[0]\n";
		print "$line[2]\n";
        } 
        else{ 
		# nuclAbundance
		print "$line[2]\n";
        } 

} 
