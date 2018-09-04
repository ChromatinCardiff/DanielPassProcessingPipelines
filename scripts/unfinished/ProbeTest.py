#!/usr/bin/python

## Written by Dr. Daniel Pass for IGEM-Cardiff 2018
## https://github.com/passdan

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
from Bio.Blast.Applications import NcbiblastnCommandline

if __name__ =='__main__':
    ## This defines the command line parameters and puts it into an array called args
    parser = argparse.ArgumentParser(description='Provide a region of DNA/RNA to be sliced to test for matches against a reference')
    parser.add_argument('-i', required=True, help='[REQUIRED] | DNA segment you want to slice for probes (FASTA FORMAT)')
    parser.add_argument('-d', required=True, help='[REQUIRED] | Blast database you want to search against (make sure it is discoverable or give full path)')
    parser.add_argument('-w', default=21, help='Window length of the probes you want to test | default = 21')
    parser.add_argument('-a', default=1, help='Number of cores you want to blast with | default = 1')
    parser.add_argument('-e', default=0.05, help='Maximum e-value for blastn-short | default = 0.05')
    parser.add_argument('-p', default=80, help='Minimum required percentage identity for blastn-short | default = 80')

    args = parser.parse_args()

def main():
    ## If the input file is a fasta (It should be!) then do this to get the main part of the name for the outputs, otherwise use the whole filename
    if args.i.endswith('.fasta'):
        outputSlug = args.i[:-6]
    elif args.i.endswith('.fa'):
        outputSlug = args.i[:-3]
    else:
        outputSlug = args.i

    probesFile =  outputSlug + "-probes.fasta"

    ## Test your blast database is correct or else fail
    testBlastDB()

    ## Open the input fasta and read it in to memory
    with open(args.i, "rU") as handle:
        for siRNA in SeqIO.parse(handle, "fasta"):
            ## Run the slicing method to get an array of the X-length fragments
            probeArray = sliceProbe(siRNA)
            ## Some nice outputs and write it to a file
            print str(len(probeArray)) + " probes were generated"

            SeqIO.write(probeArray, probesFile, "fasta")

    ## Run blast on all the probes generated
    blastOutput = blastn(probesFile, outputSlug)
    print("There were " + str(sum(1 for line in open(blastOutput))) + " matches against the blast database")
    print("\nOutputted files are:\n -- " + probesFile + "\n -- " + blastOutput)


def sliceProbe(siRNA):
    print "\n -- Reading siRNA and slicing it into windows of", args.w, "size"

    probes = []
    i = 0

    ## output the first 21 bases (or user specified) in fasta format (needed to go into blastn) then move one along base at a time, outputting each one
    while i <= len(siRNA.seq) - args.w:
        probeSeq = siRNA.seq[i:i+args.w]
        probes.append(SeqRecord(probeSeq, id=str(probeSeq), description=str(siRNA.id)))
        i += 1

    return probes

def blastn(probesFile, outputSlug):
    ## define the parameters
    print(" -- Blasting probe sequence                   " + probesFile)
    print(" -- Using BLASTN-SHORT against the database   " + args.d)
    print(" -- Requiring minimum e-value of              " + str(args.e))
    print(" -- Requiring minimum percentage identity of  " + str(args.p) )

    blast_in = (probesFile)
    blast_out = (outputSlug + "-probes.blast")

    ## Define the command. NOTE: Running blastn-short!
    # OUTPUT COLUMNS ARE: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, sequence
    blastn_cline = NcbiblastnCommandline(task="blastn-short", query=blast_in, db=args.d, perc_identity=args.p, evalue=args.e, outfmt=6, out=blast_out, num_threads=args.a)

    stdout, stderr = blastn_cline()

    blast_err_log = open("blast_err.txt", "w")
    blast_stdout_log = open("blast_stdout.txt", "w")

    blast_err_log.write(stderr)
    blast_stdout_log.write(stdout)

    return blast_out

def testBlastDB():
    try:
        open(args.d + ".nin")
    except:
        print "\n  This doesn't look like a real blast database. Did you link to the fasta file by mistake?"
        print "  If you have the fasta file but not a blast database yet, run the following command:"
        print "      makeblastdb -in myReferences.fasta -dbtype nucl -out myReferences\n"
        print "  Now you should have 4 files:"
        print "  --myReferences.fasta   (Your original input)\n  --myReferences.nhr\n  --myReferences.nin\n  --myReferences.nsq\n"
        print "  You might also have multiple of these files with numbers and a file called myReferences.nal if you've got a lot of sequences"
        print "  The script should now run with <.... -d myReferences ....>. Note that there is no suffix on that parameter"

        print "\n  ##Now quitting the script##\n"

        quit()

if __name__ == "__main__":
    main()
