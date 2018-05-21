#!/usr/bin/python

from __future__ import print_function
import os
import argparse
from Bio import SeqIO

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Import fasta file and generate x copies of each sequence based on number in header')
    parser.add_argument('-i', help='input fasta file')
    parser.add_argument('-o', help='output fasta file')
    parser.add_argument('-c', default=3, type=int, help='column number that has the multiplier value in (default is output from ampliSAS: 3)')
    parser.add_argument('-x', default=6000, type=int, help='Target number of reads (default: 5000)')
    parser.add_argument('-b', help='Artificial barcode to add to file')

    args = parser.parse_args()

def main():
    # Make an easy output file
    base=os.path.basename(args.i)
    if args.o:
        outfile = open(args.o, 'w')
    else:
        outfile = open(os.path.splitext(base)[0] + "-duplicated.fasta", 'w')

    # Do reading and multiplying
    seq_records = list(SeqIO.parse(args.i, "fasta"))

    dupFactor = findDupFactor(seq_records)

    for rec in seq_records:
        wholehead = rec.description.split('|')
        dupVal = int(wholehead[args.c].split('=')[1])
        for a in range(dupVal * dupFactor):
            print('>' + rec.id + '-' + str(a), file=outfile)
            if args.b:
                print(args.b + rec.seq, file=outfile)
            else:
                print(rec.seq, file=outfile)

def findDupFactor(seq_records):
    totalReads = 0
    for rec in seq_records:
        wholehead = rec.description.split('|')
        dupVal = int(wholehead[args.c].split('=')[1])
        totalReads += dupVal

    dupFactor = int(args.x / totalReads)

    print("Number of reads are " + str(totalReads) + " and target is " + str(args.x) + " therefore Duplication Factor is: " + str(dupFactor))

    return dupFactor

if __name__ == "__main__":
    main()
