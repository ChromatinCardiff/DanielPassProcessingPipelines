#!/usr/bin/perl
use strict;
use warnings;
use Bio::DB::Sam;

use Getopt::Std;
my %options=();
getopts('i:o:', \%options);

my $usage = "insilico_digestion_estimate.pl -i input.sam -o outfile.txt\n";


#open(OUT, '>', $options{o}) or die "$usage:$!";
print "Reading BAM file\n";

my $NATsam = Bio::DB::Sam->new(-bam  =>"$options{i}");
my @NATalignments = $NATsam->get_features_by_location();

print scalar @NATalignments . " sequences input\n";

my %nathash;
my $cut;

print "Mapping cutsites\n";
for my $a (@NATalignments){
  #print $a->qual . "\t" . $a->flag . "\n";
  if (($a->flag & 64) or ($a->flag & 80)){
    $cut = $a->start;
  #  print "$cut\n";
    ++$nathash{$a->seq_id}{$a->start};
  }elsif (($a->flag & 128) or ($a->flag & 144)){
    $cut = $a->end;
  #  print "$cut\n";
    ++$nathash{$a->seq_id}{$a->end};
  }
  # Build hash when endpoint is in bin
  #my $bin = (int($cut/ $binsize) * $binsize);
	#++$nathash{$bin};
}

print "$options{i}\n";
for my $chr (sort keys %nathash) {
    my $val = $nathash{$chr};
    my $hashsize = keys %$val;
    print "$chr is $hashsize\n";
}
print "---\n";
