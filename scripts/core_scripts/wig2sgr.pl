#!/usr/bin/perl 
use strict;
use warnings;

my $usage = "## USAGE:  wig2sgr.pl inputfile.wig outputfile.sgr\n";

my $infile = $ARGV[0];
my $outfile = $ARGV[1];

open(WIGIN, "$infile") or die $usage, $!;
open(SGROUT, ">", "$outfile") or die $usage, $!; 

my $chrom = '';
my $start = 0;
my $step = 0;
my $itt = 0;

select(SGROUT);
while(my $l = <WIGIN>){
	chomp($l);


	if ($l =~ /^fixed/){
		print STDOUT "## Procesing: $l\n";
		$itt = 0;

	        my @splitline = split(/\s+/, $l);
		
		foreach(@splitline){
			my @s = split(/=/, $_);

			if ($s[0] =~ "chrom"){
				$chrom = $s[1];
			}elsif ($s[0] =~ "start"){
				$start = $s[1];
			}elsif ($s[0] =~ "step"){
				$step = $s[1];
			}
		}
	}else{
		my $pos = $start + ($step * $itt);
		print "$chrom\t$pos\t$l\n";
		$itt++;
        } 
} 
