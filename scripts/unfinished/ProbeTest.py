#!/usr/bin/python

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
#import MOODS
import argparse
from Bio.Blast.Applications import NcbiblastnCommandline



if __name__ =='__main__':
    parser = argparse.ArgumentParser(description='Provide a region of DNA/RNA to be sliced to test for matches against a reference')
    parser.add_argument('-i', help='DNA segment you want to slice for probes (FASTA FORMAT)')
    parser.add_argument('-w', default=21, help='Window length of the probes you want to test | default = 21')
    parser.add_argument('-d', help='Blast database you want to search against (make sure it is discoverable or give full path)')
    parser.add_argument('-a', default=1, help='Number of cores you want to blast with | default = 1')

    args = parser.parse_args()


def main():
    if args.i:
        if args.i.endswith('.fasta'):
            outputSlug = args.i[:-6]
            
        with open(args.i, "rU") as handle:
            for siRNA in SeqIO.parse(handle, "fasta"):
                probeArray = sliceProbe(siRNA)
                print str(len(probeArray)) + " probes were generated"
                probesFile =  outputSlug + "-probes.fasta"
                SeqIO.write(probeArray, probesFile, "fasta")
    else:
        parser.print_help()
        quit()

    blastn(probesFile, outputSlug)


def sliceProbe(siRNA):
    print "Reading siRNA and slicing it into windows of", args.w, "size"

    probes = []
    i = 0

    while i <= len(siRNA.seq) - args.w:
        probeSeq = siRNA.seq[i:i+args.w]
        probes.append(SeqRecord(probeSeq, id=str(probeSeq), description=""))
        i += 1

    return probes

def blastn(probesFile, outputSlug):
    print("blasting probe sequence " + probesFile + " against the " + args.d + " database")

    blast_in = (probesFile)
    blast_out = (outputSlug + "-probes.blast")

    # OUTPUT COLUMNS ARE: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, sequence
    blastn_cline = NcbiblastnCommandline(task="blastn-short", query=blast_in, db=args.d, evalue=0.01, outfmt=6, out=blast_out, num_threads=args.a)

    stdout, stderr = blastn_cline()

    blast_err_log = open("blast_err.txt", "w")
    blast_stdout_log = open("blast_stdout.txt", "w")

    blast_err_log.write(stderr)
    blast_stdout_log.write(stdout)

    return blast_out


# def primer2pwm(primer):
#     """
#     Write a primer sequence as a position weight matrix.
#     """
#     # Create 4 lists of length equal to primer's length.
#     matrix = [[0] * len(primer) for i in range(4)]
#     # List of correspondance IUPAC.
#     IUPAC = {
#         "A" : ["A"],
#         "C" : ["C"],
#         "G" : ["G"],
#         "T" : ["T"],
#         "U" : ["U"],
#         "R" : ["G", "A"],
#         "Y" : ["T", "C"],
#         "K" : ["G", "T"],
#         "M" : ["A", "C"],
#         "S" : ["G", "C"],
#         "W" : ["A", "T"],
#         "B" : ["C", "G", "T"],
#         "D" : ["A", "G", "T"],
#         "H" : ["A", "C", "T"],
#         "V" : ["A", "C", "G"],
#         "N" : ["A", "C", "G", "T"]
#     }
#     # Position of nucleotides in the PWM.
#     dico = {"A" : 0,  "C" : 1, "G" : 2, "T" : 3}
#     # Read each IUPAC letter in the primer.
#     for index, letter in enumerate(primer):
#         for nuc in IUPAC.get(letter):
#             i = dico.get(nuc)
#             matrix[i][index] = 1
#     return matrix


if __name__ == "__main__":
    main()
