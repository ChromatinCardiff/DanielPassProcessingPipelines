#!/usr/bin/perl

###################################################
## Daniel Pass | github.com/passdan | March 2016 ##
###################################################

use strict;
use warnings;
use Cwd;

##############################
# Get options
##############################

use Getopt::Std;
my %options=();
getopts('i:o:', \%options);

my $usage = '## USAGE: Example:$ 3DGenomic_size_plot.pl -i infile.sam -o output.txt
				[REQUIRED]
					-i infile.sam
					-o output.txt\n';

if (!%options){
	print "~~\n$usage\n~~\n" and die;
}

#####################
# Parameter set up ##
#####################

my $infile = $options{i};
my $outfile = $options{o};




#################
# Main Program ##
#################

my %sizehash;
#%covhash;

open(IN, "$infile") || die "Unable to open $infile: $!";

while(<IN>){
	my @line = split('\t',$_);

	if ($line[8] > 0){
		my $posinc = $line[3];

		while ($posinc <= $line[7]){
			if (exists $sizehash{$posinc}){
				$sizehash{$posinc} .= ",$line[8]";
			}else{
				$sizehash{$posinc} = "$line[8]";
			}
			#$covhash{$posinc} += 1;
			$posinc++;
		}
	}
}

open(OUT, '>', $outfile) or die "Unable to open $outfile: $!";
select OUT;
# Export bin counts per chromosome;
print "Pos\tSizes\n";
for my $baseN (sort {$a <=> $b} keys %sizehash){
	my @cov = split(',',$sizehash{$baseN});
	foreach (@cov){
		my $round = int($_ / 10) *10;
		print "$baseN\t$round\n";
	}
    #print OUT "$baseN\t" . scalar @cov ."\t$sizehash{$baseN}\n";
	#print OUT "$baseN\t$sizehash{$baseN}\n";
}
