#!/usr/bin/perl

# Daniel Pass | daniel.antony.pass@googlemail.com | 16/04/2015

use strict;
use warnings;
use Bio::DB::Sam;

##############################
# Get options
##############################
use Getopt::Std;
my %options=();
getopts('i:o:m:t:b:IRNh', \%options);
my $sep = "~~~\n";

my $usage = "\n\
AbsenceMapper is used to identify regions of a mapped genome which don't demonstrate any endpoints of sequencing reads.\ This could be useful when determining differential degradation from MNase-seq experiments, but also identifying unmappable regions or missed sections in low coverage sequencing.\n\
USAGE: absencemapper.pl -i infile.sam [REQUIRED]\n\
[OPTIONAL]\n \
 -o [outdir | default: .]\n \
 -m [min region size to report | Default: 100]\n \
 -t [Threshold for reads bound to cut site to count (options: LQ/MEDIAN/UQ or integer) | default: LQ]\n \
 -b [binsize | default: 5]\n  -n [Noise level to ignore at gap determination and mag calculation | default 1]\n \
 -I [Ignore singletons (one count bins) when calculating thresholds]\n\
 -R [Suppres SAM parsing/cutfile generation]\n \
 -N [Do not normalise cut intensity by dataset mean | default: False]\n\n \
 -h Prints this helpful parameter list and exits\n\n\
NOTE: [-R] reduces processing time but assumes files are in correct outdir and no change to names\n";

if ((!%options) or ($options{h})){
	print "$usage$sep" and die;
}

# Infile name
my $infile = $options{i};
(my $prefix = $infile) =~ s/\.sam//;

# Assign defaults
my ($outdir, $cutoutfile, $gapoutfile, $mingap, $thresh, $binsize, $MININT);
my @thresholds;

# Out directory
if(exists $options{o}){
	$outdir = $options{o};
	mkdir $outdir;
}else{
	$outdir = ".";
}

# Min gap size to report
if(exists $options{m}){
	$mingap = $options{m};
}else{
	$mingap = 100;
}

# Min magnitude of cut size to count
if(exists $options{t}){
	@thresholds = split(',', $options{t});
}else{
	push(@thresholds, "LQ");
}

# Bin size for cutsites
if(exists $options{b}){
	$binsize = $options{b};
}else{
	$binsize = 5;
}

## RUN PROGRAM

# Print parameters to screen
print $sep;
print "Input____________________________________________$infile\n";
print "Output directory_________________________________$outdir\n";
print "Minimum distance between cuts to report__________$mingap\n";
print "Minimum cuts at site required____________________@thresholds\n";
print "Size of bin for each cutsite_____________________$binsize\n";
print "Ignore singletons?_______________________________"; if(exists($options{I})){print "Yes\n";}else{print "No\n";};
print "Normalise data by readcount per chromosome?______"; if(!exists($options{N})){print "Yes\n";}else{print "No\n";};
print "Skip SAM parsing (Jump straight to gap calcs)?___"; if(exists($options{R})){print "Yes\n";}else{print "No\n";};
print $sep;

# Main hash file!
my %cuthash;
my $cuthashsize;
my $desinghashsize;

## Parse SAM file

if(!exists $options{R}){
	print "Reading SAM file\n$sep";
  my $sam = Bio::DB::Sam->new(-bam  =>"$infile");

  while (my $align = $sam->
  my $iterator     = $segment->features(-iterator=>1);
  while (my $align = $iterator->next_seq) { do something }
	# open(IN, "<", "$infile") || die "Unable to open $infile: $!";
  #
	# my (@metadata, $sum);
  #
  # print "Filtering non-mapped reads\n";
  #
	# while(my $line = <IN>){
	# 	# Ignore metadata
	# 	if($line =~ /^@/){
	# 		push(@metadata, $line);
	# 	}else{
	# 		# split line by delimiter and store elements in an array
	# 		my @line = split('\t',$line);
	# 		# Ignore '*' failures
	# 		if ($line[2] ne '*'){
	# 			my $bin = (int($line[7]/ $binsize) * $binsize);
	# 			# Build hash when endpoint is in bin
	# 			++$cuthash{$line[2]}{$bin};
	# 		}
	# 	}
	# }
	close(IN);

	## Full hashsize
	for my $chr (sort keys %cuthash){
		my $val = $cuthash{$chr};
		$cuthashsize += keys %$val;
	}

	## Normalise reads by chromosome (if -N passed)

	if(!exists $options{N}){
		for my $chr (sort keys %cuthash){
			print "Normalising $chr by number of reads ($cuthashsize)\n";
			for my $cutsite (keys %{$cuthash{$chr}}){
				$cuthash{$chr}{$cutsite} = ($cuthash{$chr}{$cutsite} / ($cuthashsize / 1000000));
			}
		}
		print $sep;
	}

	## Ignoring singletons?
	if(exists $options{I}){
		print "Dropping singletons (bins with only 1 cut)\n";
		my %desinghash;
		my $decount;

		for my $chr (keys %cuthash){
			while (my ($pos, $value) = each %{ $cuthash{$chr} }){
				#	print "$key1\t$key2\t$value\n";
				if ($value > 1){
					$desinghash{$chr}{$pos} = $value;
					$decount++;
				}
			}
		}
		%cuthash = %desinghash;

		print "Full dataset cutsite count:	$cuthashsize\n";
		print "After singleton removal:	$decount\n$sep";
	}

	## Generate cutsite file
  my %FH;

	foreach my $val (@thresholds){
		$cutoutfile = "$prefix" . "_cuts_$val.txt";
		print "opening $cutoutfile\n";
		open(my $threshhandle, ">", "$outdir/$cutoutfile") or die "Unable to open $cutoutfile: $!";
		$FH{$val} = $threshhandle;
	}
	print $sep;

	# Export bin counts per chromosome;
	print "Frequency of cuts per site, Quartile range (x10^-6):\n";

	for my $chr (sort keys %cuthash){
    my %TH;

		# calculate quartile limit per chromosome
		my @sortvals = sort {$a <=> $b} values %{$cuthash{$chr}};
		my $max = scalar @sortvals;

    $TH{"LQ"} = $sortvals[int(0.25 * $max)];
		$TH{"MED"} = $sortvals[int(0.5 * $max)];
		$TH{"UQ"} = $sortvals[int(0.75 * $max)];

		foreach my $val (@thresholds){
			if ($val =~ /^\d+$/){
				$TH{$val} = $sortvals[int(($val/100) * $max)];
			}
		}

		my @keys = sort { $TH{$a} <=> $TH{$b} } keys(%TH);
		my @vals = @TH{@keys};

		printf <%-12s>, $chr;
		foreach my $key (@keys) {
		    print "\t$key:\t";
				printf ("%.3f", $TH{$key});
				print "\t|";
		}
		print "\n";

		for my $cutsite (sort {$a <=> $b} keys %{$cuthash{$chr}}){
			foreach my $val (@thresholds){
#				if ($val =~ /^\d+$/){
#					$MININT = ###;
#					print $MININT;
#				}else{
					$MININT = $TH{$val};
#				}

				if ($cuthash{$chr}{$cutsite} >= $MININT){
					$cuthash{$val}{$chr}{$cutsite} = $cuthash{$chr}{$cutsite};
					print {$FH{$val}} "$chr\t$cutsite\t$cuthash{$chr}{$cutsite}\n";
				}
			}
		}
	}
	print STDOUT $sep;
}
select STDOUT;

## Regenerate cuthash if skipping the SAM parsing step

if(exists $options{R}){
	print "Skipping SAM parsing, regenerating hash from cuts file\n";
	foreach my $val (@thresholds){
		print "opening $prefix" . "_cuts_$val.txt\n";
		open(CUTSIN, "<", "$outdir/$prefix" . "_cuts_$val.txt") || die "Unable to open $outdir/$prefix" . "_cuts_$val.txt: $!";
		while (my $line = <CUTSIN>){
			my @line = split('\t',$line);
			$cuthash{$val}{$line[0]}{$line[1]} = $line[2];
		}
		close(CUTSIN);
	}
	print $sep;
}
#############################
## Identify gapped regions ##
#############################
my %GH;

foreach my $val (@thresholds){
	$gapoutfile = "$prefix" . "_gaps_$val.gff";
	print STDOUT "opening $gapoutfile\n";
	open(GAPOUT, ">", "$outdir/$gapoutfile") or die "Unable to open $gapoutfile: $!";
  select(GAPOUT);
	my ($start,$end) = 0;
	my $gapcounttotal;

	print "##gff-version 3\n";
	print STDOUT "Number of cut-free regions at $val threshold:\n";
	for my $chr (sort keys %{$cuthash{$val}}){
		print STDOUT "--$chr:\t";
		my $gapcountchr;
		for my $cutsite (sort {$a <=> $b} keys %{$cuthash{$val}{$chr}}){
			$gapcounttotal++;
			$end = $cutsite;
			my $gap = $end - $start;

			if ($gap > $mingap){
				print "$chr\tProtDNA\texon\t$start\t$end\t.\t.\t.\n";
				$gapcountchr++;
			}
			$start = $cutsite;
		}
		print STDOUT "$gapcountchr\n";
	}
  close(GAPOUT);
	print STDOUT "\n$infile reports $gapcounttotal cut free regions at $val threshold\n";
	print STDOUT "Output:\t$gapoutfile\n$sep";
}

select STDOUT;
print "Job done!\n";
