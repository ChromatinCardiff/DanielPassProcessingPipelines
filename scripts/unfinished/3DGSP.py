#!/usr/bin/python

###################################################
## Daniel Pass | github.com/passdan | March 2016 ##
###################################################

import argparse
#from Bio import SeqIO
import pysam


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Import 3DG files for plotting')
    parser.add_argument('-i', required=True, help='Input Bam file')
    parser.add_argument('-o', help='Output prefix')
    parser.add_argument('-b', help='Regions to extract (Bed file)')
    parser.add_argument('-B', help='Individual set of co-ordinates in Chr1,1000,1500 format (include commas)')
    parser.add_argument('-e', default=0, type=int, help='Integer for upstream/downstream extension')

    args = parser.parse_args()

def main():
    samfile = pysam.AlignmentFile(args.i, "rb")

    for read in samfile.fetch('Chr1', 100000, 100800):
        #line = read.split('\t')
        print read.query_alignment_start(), read.query_alignment_end(), read.query_length()

    samfile.close()

if __name__ == "__main__":
    main()
