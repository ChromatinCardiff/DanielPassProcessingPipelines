#!/usr/bin/python

#################################################
## Daniel Pass | github.com/passdan | Oct 2016 ##
#################################################

# Basic processing
import argparse
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='merge multiple wigfiles')
    parser.add_argument('-i', help='Inputfiles (comma seperated)')
    parser.add_argument('-o', help='outputfile')
    parser.add_argument('-m', default="mean", help='method to merge (default: mean)')
    args = parser.parse_args()

def main():
    sampleList = []
    bCount = 0
    # If merging multiple bamfiles
    if args.i:
        wiglist = args.i.split(',')
        for file in wiglist:
            #print "Reading in file:", file
            with open(file, 'r') as infile:
                f = infile.read().splitlines()
                sampleList.append(f)
    else:
        print "Need an input file! or use -h to see options"
        quit()

    outfile = open(args.o, 'w')

    sumwig = []
    fileinc = 0

    for file in sampleList:
        lineinc = 0
        for line in file:
            if re.match('fix', line):
                if fileinc == 0:
                    sumwig.append(line)
            else:
                if line.strip():
                    value = float(line)
                if fileinc == 0:
                    sumwig.append(value)
                else:
                    sumwig[lineinc] = sumwig[lineinc] + value
            lineinc += 1

        fileinc += 1

    for line in sumwig:
        if re.match('fix', str(line)):
            outfile.write(''.join(line + "\n"))
        else:
            outfile.write(''.join((str(float(line) / float(fileinc))+"\n")))

def updateHeadLine(line):
    splitLine = line.split(' ')
    inputInterval = splitLine[3].split('=')[1]
    newInterval = int(args.n) * int(inputInterval)
    newLine = (' '.join((splitLine[0],splitLine[1],splitLine[2], 'step=' + str(newInterval), 'span=' + str(newInterval), '\n')))
    return newLine

if __name__ == "__main__":
    main()
