#!/usr/bin/python


#################################################
## Daniel Pass | github.com/passdan | Nov 2016 ##
#################################################

# Basic processing
import argparse
import re
import operator


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='count number of intervals in genome')
    parser.add_argument('-i', help='inputfile')
    parser.add_argument('-o', help='outputfile')
    parser.add_argument('-n', default=10000,help='size of interval (bp)')
    args = parser.parse_args()

def main():
    with open(args.i, 'r') as infile:
        genes = infile.read().splitlines()
    outfile = open(args.o, 'w')

    chromosomes = { 'Chr1':30427680,
                    'Chr2':19698290,
                    'Chr3':23459840,
                    'Chr4':18585060,
                    'Chr5':26975510,
                    }

    intSize = int(args.n)
    intCounter = 0
    intStart = 0

    chromoDB = buildChromoDB(chromosomes,intSize)

    for line in genes:
        #print line
        chrom = line.split('\t')[0]
        gStart = line.split('\t')[3]
        genebin = ((int(int(gStart)/intSize))*intSize + 1)
        #print chrom, str(gStart), genebin
        chromoDB[(chrom,genebin)] += 1
        #print chrom, genebin, chromoDB[(chrom,genebin)]

    printGenomeDB(chromoDB, outfile)
    #print chromoDB


def buildChromoDB(chromosomes, intSize):
    i = 1
    chromoDB ={}
    for key in chromosomes:
        while i < chromosomes[key]:
            chromoDB[(key,i)] = 0
            i += intSize
        i = 1

    #printGenomeDB(chromoDB)
    return chromoDB

def printGenomeDB(DB, outfile):
    for coords, value in sorted(DB.items(), key=operator.itemgetter(0)):
        outfile.write(''.join((coords[0] + "\t" + str(coords[1]) + "\t" + str(value) + "\n")))

if __name__ == "__main__":
    main()
