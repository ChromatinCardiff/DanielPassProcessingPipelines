#!/usr/bin/python

##################################################
## Daniel Pass | github.com/passdan | June 2016 ##
##################################################

import argparse
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
import re

description = """---------------
Using info from: http://www.clontech.com/GB/Products/Genome_Editing/CRISPR_Cas9/Resources/Designing_sgRNA
Use the following guidelines to select a genomic DNA region that corresponds to the crRNA sequence of the sgRNA:
- The 3' end of the DNA target sequence must have a proto-spacer adjacent motif (PAM) sequence (5'-NGG-3').
- The 20 nucleotides upstream of the PAM sequence will be your targeting sequence (crRNA) and Cas9 nuclease will cleave approximately
3 bases upstream of the PAM.
- The PAM sequence itself is absolutely required for cleavage, but it is NOT part of the sgRNA sequence and therefore should not
be included in the sgRNA.
- The target sequence can be on either DNA strand.

Tips for designing sgRNAs: Through our own experience, we have identified additional tips for designing sgRNAs.
We have found that the best sgRNAs for several tested genes have a G at position 1 and an A or T at position 17.
---------------"""


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='test')
    parser.add_argument('-i', help='Input test')
    parser.add_argument('-o', default="cas9-output", help='Output handle')
    parser.add_argument('-B', action='store_true', help='Pass argument to Blast the results')
    parser.add_argument('-d', default='ecoli-k12-genome', help='Database to Blast the results against')

    args = parser.parse_args()

def main():
    print description
    outfasta = open(args.o + ".fasta", 'w')
    outtable = open(args.o + "-table.txt", 'w')


    for seq_record in SeqIO.parse(args.i, "fasta"):
        i = 0
        targets =[]
        targetDict = {}
        fwdDict = {}
        revDict = {}

        print "Building 20bp dictionary."
        while i < len(seq_record):
            targets.append(seq_record.seq[i:(i+23)])
            targetDict[i] = seq_record.seq[i:(i+23)]
            i += 1

        print "Testing if regions match cas9 design guidlines."
        outtable.write(' '.join(("Gene:", seq_record.description, "\n")))
        outtable.write("Pos  : 5' - 20bp sgRNA      | PAM |  Spacer | PAM | RevComp 5'-20bp sgRNA\n")
        outtable.write("========================================================================\n")

        for pos in sorted(targetDict.keys()):
            t = targetDict[pos]
            if re.match('^G.{15}[A|T].{3}.GG$', str(t)) is not None:
                if re.search('N', str(t)) is not None:
                    print '[NB: ambiguous bases present]'
                else:
                    a = t[0:20], "|", t[20:23]
                    fwdDict[pos] = ' '.join(str(i) for i in a)

        for pos in sorted(targetDict.keys()):
            t = targetDict[pos].complement()
            if re.match('^GG..{3}[A|T].{15}G', str(t)) is not None:
                if re.search('N', str(t)) is not None:
                    print '[NB: ambiguous bases present]'
                else:
                    a = t[0:3], "|", t[3:23]
                    revDict[pos] = ' '.join(str(i) for i in a)

        matchInt = 0

        for fwdpos in sorted(fwdDict.keys()):
            for revpos in sorted(revDict.keys()):
                spacer = (revpos + 3) - (fwdpos + 20)
                if (spacer < 100) and (spacer > 5):
                    outtable.write(' '.join((str(fwdpos).zfill(4),":", fwdDict[fwdpos], "--", str(spacer), "bp--", revDict[revpos], "\n")))
                    matchInt += 1
                    outfasta.write(''.join((">", str(i + matchInt) + "-forward | FwdPos:", str(fwdpos).zfill(4), "\n")))
                    outfasta.write(''.join((fwdDict[fwdpos][0:20], "\n")))
                    outfasta.write(''.join((">", str(i + matchInt) + "-reverse | FwdPos:", str(fwdpos).zfill(4), "\n")))
                    outfasta.write(''.join((revDict[revpos][6:26], "\n")))

    if args.B is True:
        print "Running blastn for all generated sequences."
        blastn_cline = NcbiblastnCommandline(query=args.o + ".fasta", db=args.d, outfmt=7, out=args.o + ".bln")
        stdout, stderr = blastn_cline()

    print "All done, see files:", args.o + "-table.txt,", args.o + ".fasta and", args.o + ".bln"

if __name__ == "__main__":
    main()
