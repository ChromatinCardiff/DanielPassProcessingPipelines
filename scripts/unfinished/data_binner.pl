#!/usr/bin/perl
use strict;
use warnings;

##############################
# Get options
##############################
use Getopt::Std;
my %options=();
getopts('i:o:b:m:', \%options);

my $usage = '## USAGE: Example:$ data_binner.pl -i infile.txt -o outfile.sgr
		[REQUIRED]
			-i infile.txt
			-o outdirectory
		[OPTIONAL]
			-b INT | binsize (default:10)
			-m counting-method | options: SUM or AVG (default: AVG)\n';

if (!%options){
	print "~~\n$usage\n~~\n" and die;
}

################################################################################
# Parameter set up
################################################################################

my $infile = $options{i};
my $outfile = $options{o};
my ($binwidth, $method);

if(exists $options{b}){
	$binwidth = $options{b};
}else{
	$binwidth = 10;
}

if(exists $options{m}){
	$method = $options{m};
}else{
	$method = "AVG";
}

################################################################################
# main program
###############################################################################

open(IN, "$infile") or die "Unable to open $infile: $!";
open (OUTFILE, '>', $outfile) or die "Unable to open outfile: $!";

my %bintotal;
my $i;

while(<IN>){
	chomp($_);
	# split line by delimiter and store elements in an array
	my @line = split('\t',$_);

	my $bin = (int($line[1]/$binwidth))*$binwidth;
	#print "Chr: $line[0]\tbin:$bin\tPos:$line[3]\n";

	if (exists $bintotal{$line[0]}{$bin}){
		$bintotal{$line[0]}{$bin} += $line[2];
	}else{
		$bintotal{$line[0]}{$bin} = $line[2];
	}

	$i++;
	if($i % 1000000 == 0){
		print "$i lines processed\n";
	}
}
close(IN);

# Export bin counts per chromosome;
for my $chr (sort keys %bintotal){
	if($method eq "SUM"){
		print "Processing $chr\n";
		for my $binout (sort {$a <=> $b} keys %{$bintotal{$chr}}){
			print OUTFILE "$chr\t$binout\t$bintotal{$chr}{$binout}\n";
		}
	}elsif($method eq "AVG"){
		print "Processing $chr\n";
		for my $binout (sort {$a <=> $b} keys %{$bintotal{$chr}}){
			my $value = $bintotal{$chr}{$binout} / 10;
			print OUTFILE "$chr\t$binout\t$value\n";
		}
	}
}
