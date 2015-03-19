#!/usr/bin/perl
use strict;
use warnings;

my $usage = "## USAGE:  filter_csv.pl infile.txt [column] [min value]\n";

my $infile = $ARGV[0] or die $usage;
my $column = $ARGV[1] or die $usage;
my $filter = $ARGV[2] or die $usage;

open INFILE, "<", "$infile" or die "$usage | check $infile\n";

my $n=0;

while (my $line = <INFILE>){
	$n++;
	if ($n == 1){
		print $line;
	}else{
		my @l = split(/\t/, $line);
		if ($l[$column] > $filter){
			print $n . "_" . $line;
		}
	}
}

