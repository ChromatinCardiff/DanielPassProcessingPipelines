#!/usr/bin/perl
use strict;
use warnings;

my $usage = "## USAGE:  feature_mapping.pl chromatin_map.sgr feature_map.data [START pos column] [END POS column] new_map.txt\n";

my $map = $ARGV[0] or die $usage;
my $features = $ARGV[1] or die $usage;
my $start = $ARGV[2] or die $usage;
my $end = $ARGV[3] or die $usage;
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
		if ($l[1] > $f[$start - 1] && $l[1] < $f[$end - 1]){
			if ($flag == 0){
				$line .= "\t$f[1]";
				$flag = 1;
			}else{
				$line .= ",$f[1]";
			}
		}
	}
	$flag = 0;
	print OUTDIR "$line\n";
}
