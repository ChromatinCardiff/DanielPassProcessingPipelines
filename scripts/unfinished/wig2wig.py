#!/usr/bin/python


#################################################
## Daniel Pass | github.com/passdan | Oct 2016 ##
#################################################

# Basic processing
import argparse


if __name__ = '__main__':
    parser = argparse.ArgumentParser(description='change interval size for wigfile')
    parser.add_argument('-i', help='inputfile')
    parser.add_argument('-o', help='outputfile')
    parser.add_argument('-n', help='new interval size in basepairs')

    args = parser.parse_args()

def main():
    inflie = open(args.i)
    




if __name__ == "__main__":
    main()
