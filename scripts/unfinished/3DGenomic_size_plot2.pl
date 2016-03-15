#!/usr/bin/perl

###################################################
## Daniel Pass | github.com/passdan | March 2016 ##
###################################################

use strict;
use warnings;
use Bio::DB::Sam;

##############################
# Get options
##############################

use Getopt::Std;
my %options=();
getopts('i:o:b:B:e:', \%options);

my $usage = '## USAGE: Example:$ 3DGenomic_size_plot.pl -i infile.sam -o output.txt
				[REQUIRED]
					-i infile.sam
					-o output.txt
				[OPTIONAL]
					-b regions.bed {read bed file with regions to extract}
					-B Chr1,100000,100500 {Individual co-ordinates for extraction, commas required}
					-e INTEGER {bases up/downstream to also inculde in region extract}\n';

if (!%options){
	print "~~\n$usage\n~~\n" and die;
}

#####################
# Parameter set up ##
#####################

my $infile = $options{i};
my $outfile = $options{o};

my @regions;
if (exists = $options{B}){
	@regions = $options{B};
}elsif (exists = $options{b}){
	open(BED, $options{b}) || die "Unable to open " . $options{b} . ": $!";
	while(<BED>){
		push(@regions, $_);
	}
}


my @

my $sam = Bio::DB::Sam->new(-bam  => $infile) || die "Unable to open $infile: $!";
my @alignments = $sam->get_features_by_location(-start=>)


####################
# Extract regions ##
####################





#################
# Main Routine ##
#################

sub REGEXTRACT {
	my %sizehash;
	#%covhash;

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
}
