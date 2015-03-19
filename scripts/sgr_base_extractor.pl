#!/usr/bin/perl
use strict;
use warnings;

my $usage = "## USAGE:  sgr_base_extractor.pl chromatin_map.sgr feature_map.csv [upstream amount] [downstream amount] new_map.txt\n";

my $map = $ARGV[0] or die $usage;
my $features = $ARGV[1] or die $usage;
my $upstream = $ARGV[2] or die $usage;
my $downstream = $ARGV[3] or die $usage;
my $newmap = $ARGV[4] or die $usage;

open CMTN, "<", "$map" or die "$usage | check $map\n";
my @cmtn = (<CMTN>);
open FEAT, "<", "$features" or die "$usage | check $features\n";
my @feat = (<FEAT>);
open OUTDIR, ">", "$newmap" or die "$usage | check $newmap\n";

foreach my $line (@cmtn){
	chomp($line);
	my $flag = 0;

	foreach my $feat (@feat){
		my @l = split(/\t/, $line);
		my @f = split(/\t/, $feat);

#		print "is $l[1]	between $f[2] and $f[3]\n?";
		if ($l[1] > ($f[5] - $downstream) && $l[1] < ($f[5] + $upstream)){
			print OUTDIR "$line\t" . $l[1] - $f[5] . "$f[1]\t$f[6]\t$f[7]\t$f[8]\t$f[10]\n";
			
		}
	}
	$flag = 0;

}
