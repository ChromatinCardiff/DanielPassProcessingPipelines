#!/usr/bin/python


#################################################
## Daniel Pass | github.com/passdan | Oct 2016 ##
#################################################

# Basic processing
import argparse
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='change interval size for wigfile')
    parser.add_argument('-i', help='inputfile')
    parser.add_argument('-o', help='outputfile')
    parser.add_argument('-n', help='number of intervals to merge together (average by mean)')
    args = parser.parse_args()

def main():
    with open(args.i, 'r') as infile:
        wig = infile.read().splitlines()
    outfile = open(args.o, 'w')

    mergeInterval = int(args.n)
    averageCounter = 0
    inc = 0

    for line in wig:
        if re.match('fix', line):
            if inc > 0:
                outfile.write(str(averageCounter / int(inc)) + '\n')
            newline = updateHeadLine(line)
            outfile.write(str(newline))
            inc = 1
            averageCounter = 0
        else:
            if inc <= mergeInterval:
                averageCounter += float(line)
                inc += 1
            else:
                outfile.write(''.join((str(averageCounter / int(mergeInterval)) + '\n')))
                inc = 0
                averageCounter = 0

def updateHeadLine(line):
    splitLine = re.split(' |=', line)
    newInterval = int(args.n) * int(splitLine[6])
    splitLine[6] = newInterval
    splitLine[8] = newInterval
    jLine = ' '.join(map(str, splitLine)) + "\n"
    kLine = jLine.replace(r"chrom ", "chrom=")
    lLine = kLine.replace(r" 1", "=1")
    return lLine


if __name__ == "__main__":
    main()
