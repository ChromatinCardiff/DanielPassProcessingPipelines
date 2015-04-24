#!/usr/bin/perl

# Daniel Pass daniel.antony.pass@googlemail.com 16/04/2015

use strict;
use warnings;

##############################
# Get options
##############################
use Getopt::Std;
my %options=();
getopts('i:o:m:t:b:IRN', \%options);
my $sep = "~~~\n";

my $usage = "\nAbsenceMapper is used to identify regions of a mapped genome which don't demonstrate any endpoints of sequencing reads. This could be useful when determining differential degradation from MNase-seq experiments, but also identifying unmappable regions or missed sections in low coverage sequencing.\n\
USAGE: absencemapper.pl -i infile.sam [REQUIRED]\n[OPTIONAL]\n  -o [outdir | default: .]\n  -m [min gap to report | Default: 50]\n  -t [Threshold for cut site to count (options: LQ/MEDIAN/UQ or integer) | default: LQ]\n  -b [binsize | default: 5]\n  -n [Noise level to ignoreIgnore singletons at gap determination and mag calculation]\n  -R [Suppres SAM parsing/cutfile generation]\n  -N [Do not normalise cut intensity by dataset mean (default: Normalise)]\n\n  -h Prints this helpful usage list and exits\n\nNOTE: [-R] reduces processing time but assumes files are in correct outdir and no change to names\n";

if (!%options){
	print "$usage$sep" and die;
}

if ($options{h}){
	print "$usage$sep" and die;
}

# Infile name
my $infile = $options{i};
(my $prefix = $infile) =~ s/\.sam//;

# Assign defaults
my ($outdir, $cutoutfile, $gapoutfile, $mingap, $minmag, $binsize);

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
	$mingap = 50;
}

# Min magnitude of cut size to count
if(exists $options{t}){
	$minmag = $options{t};
}else{
	$minmag = "LQ";
}

# Bin size for cutsites
if(exists $options{b}){
	$binsize = $options{b};
}else{
	$binsize = 5;
}

# Print parameters to screen
print $sep;
print "Input____________________________________________$infile\n";
print "Output directory_________________________________$outdir\n";
print "Minimum distance between cuts to report__________$mingap\n";
print "Minimum cuts at site required____________________$minmag\n";
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
	print "Parsing SAM file\n$sep";
	open(IN, "<", "$infile") || die "Unable to open $infile: $!";

	my (@metadata, $sum);

	while(my $line = <IN>){
		# Ignore metadata
		if($line =~ /^@/){
			push(@metadata, $line);
		}else{
			# split line by delimiter and store elements in an array
			my @line = split('\t',$line);
			# Ignore '*' failures
			if ($line[2] ne '*'){
				my $bin = (int($line[7]/ $binsize) * $binsize);
				# Build hash when endpoint is in bin
				++$cuthash{$line[2]}{$bin};
			}
		}
	}
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
			while (my ($pos, $value) = each %cuthash{$chr}){
				#	print "$key1\t$key2\t$value\n";
				if ($value > 1){
					$desinghash{$chr}{$pos} = $value;
					$decount++;
				}
			}
		}
		%cuthash = %desinghash;

		print "Full dataset:  $cuthashsize\n";
		print "Post removal:  $decount\n$sep";
	}

	## Generate cutsite file

	my ($LQ,$MED,$UQ,$MININT);

	if ($minmag =~ /^d+$/){
		$MININT = $minmag;
	}

	$cutoutfile = "$prefix" . "cuts_$minmag.txt";
	open(CUTOUT, '>', "$outdir/$cutoutfile") or die "Unable to open $cutoutfile: $!";
	select(CUTOUT);
	# Export bin counts per chromosome;
	for my $chr (sort keys %cuthash){

		# calculate quartile limit per chromosome
		my @sortvals = sort {$a <=> $b} values %{$cuthash{$chr}};
		my $max = scalar @sortvals;

		#print STDOUT "max $max\nLQ: " . int(0.25 * $max) . "\nsortvals[" . int(0.25 * $max) . "]: " . $sortvals[32094] . "\n";
		my $LQ = $sortvals[int(0.25 * $max)];
		my $MED = $sortvals[int(0.5 * $max)];
		my $UQ = $sortvals[int(0.75 * $max)];

		if ($minmag =~ /LQ/){
			$MININT = $LQ;
		}elsif ($minmag =~ /MED/){
			$MININT = $MED;
		}elsif ($minmag =~ /UQ/){
			$MININT = $UQ;
		}
		print STDOUT "$chr\tQuartile Range (x10^-6):\tLQ:" . sprintf("%.3f",$LQ) . "\tMEDIAN:" . sprintf("%.3f",$MED) . "\tUQ:" . sprintf("%.3f",$UQ) . "\tSelected: $minmag(" . sprintf("%.3f",$MININT) . ")\n";

	  for my $cutsite (sort {$a <=> $b} keys %{$cuthash{$chr}}){
			if ($cuthash{$chr}{$cutsite} >= $MININT){
	      print "$chr\t$cutsite\t$cuthash{$chr}{$cutsite}\n";
			}
	  }
	}
	close(CUTOUT);
	print STDOUT $sep;
}
select STDOUT;

## Regenerate cuthash if skipping the SAM parsing step

if(exists $options{R}){
	print "Skipping SAM parsing, regenerating hash from cuts file\n";
	open(CUTSIN, "<", "$outdir/$cutoutfile") || die "Unable to open $outdir/$cutoutfile: $!";
  while (my $line = <CUTSIN>){
    my @line = split('\t',$line);
		$cuthash{$line[0]}{$line[1]} = $line[2];
	}
	print $sep;
}

$gapoutfile = "$prefix" . "_gaps_$minmag.txt";
## Identify gapped regions
open(GAPOUT, '>', "$outdir/$gapoutfile") or die "Unable to open $outdir/$gapoutfile: $!";

my ($previous,$current) = 0;

select(GAPOUT);
my $gapcounttotal;
print "Chr\tStart\tEnd\tSize\n";
for my $chr (sort keys %cuthash){
	print STDOUT "Processing $chr for cut-free regions\n";
	my $gapcountchr;
	for my $cutsite (sort {$a <=> $b} keys %{$cuthash{$chr}}){
		$gapcounttotal++;
		$current = $cutsite;
	  my $gap = $current - $previous;

		if ($gap > $mingap){
	    print "$chr\t$previous\t$current\t$gap\n";
			$gapcountchr++;
	  }
	  $previous = $cutsite;
	}
	print STDOUT "$chr\tCut free regions: $gapcountchr\n";
}
select STDOUT;
print $sep;
print "$infile reports $gapcounttotal cut free regions\n";
print "Output:\nCutsite file: $cutoutfile\nIdentified gaps file: $gapoutfile\n";
print $sep;
print "Job done!\n";
