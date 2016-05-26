#!/usr/bin/python

#################################################
## Daniel Pass | github.com/passdan | May 2016 ##
#################################################

import argparse
import re


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter gtf to remove annotations from different chromosomes')
    parser.add_argument('-i', help='input gtf)')

    args = parser.parse_args()


def main():

    infile = open(args.i, "r")
    head = infile.readline()

    trans = 0
    chrom = 0
    for f in infile:
        fsplit = f.split()

        # num = 1
        # for i in fsplit:
        #     print num, i
        #     num += 1
        #print trans, fsplit[9]
        #print chrom, fsplit[0]

        if chrom == fsplit[0] and trans == fsplit[9]:
            print f.rstrip()
            chrom = fsplit[0]
        if trans != fsplit[9]:
            print f.rstrip()
            chrom = fsplit[0]


#        chrom = fsplit[0]
        trans = fsplit[9]


if __name__ == "__main__":
    main()
