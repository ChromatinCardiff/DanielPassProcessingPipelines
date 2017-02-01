#!/usr/bin/python

#################################################
## Daniel Pass | github.com/passdan | Nov 2016 ##
#################################################

# Still under development! Keep an eye out for updates! #

# Basic processing
import argparse
import pysam
import re
from collections import defaultdict
from datetime import datetime
startTime = datetime.now()
# Numbers
import numpy as np
from scipy import stats
# Graphing
import plotly
import plotly.graph_objs as go

## Define some parameters IF YOU WANT TO SAVE TYPING! ##
SAMdir = "/home/sbi6dap/RawData/ALD/MNase-seq/BAMs/"
GTFFile = "/home/sbi6dap/Projects/REFDB/Araport11_GFF3_genes_transposons.201606.gtf"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Import Bam file and extract defined region for 2D or 3D plotting')
    parser.add_argument('-i', help='Input a Bam file or a list of comma-separated Bam files to merge (format: one.bam,two.bam,three.bam)')
    parser.add_argument('-d', default=SAMdir, help='Optional: Directory where bamfiles are (saves on typing!)')
    parser.add_argument('-o', help='Output prefix')
    parser.add_argument('-B', help='One individual set of co-ordinates in bed format (Chr1:1000-1510)')
    parser.add_argument('-b', default=False, help='Bed file of regions to extract (Bed file)')
    parser.add_argument('-g', help='Gene transcript to get positional information from | REQUIRES GTF FILE (-G)')
    parser.add_argument('-G', default=GTFFile, help='gtf file to get gene model information from | REQUIRES GENE/TRANSCRIPT ID')
    parser.add_argument('-e', default=500, type=int, help='Integer for upstream/downstream extension')
    parser.add_argument('-s', default=50, type=int, help='Minimum size of particle')
    parser.add_argument('-S', default=1000, type=int, help='Maximum size of particle')
    parser.add_argument('-q', default=10, type=int, help='Integer for rounding counts (fragment size, default: 10)')
    parser.add_argument('-r', default=10, type=int, help='Integer for rounding counts (bases, default: 50)')
    parser.add_argument('-p', default="2D", help='Choose plot type (2D or 3D, default:2D)')
    parser.add_argument('-z', default=1000, help='Change zmax for plot (Mainly to look at low abundance things) (Default 1000, \'None\' allows dynamic assignment)')
    parser.add_argument('-t', action='store_true', help='Use default co-ordinates for testing')

    args = parser.parse_args()

def main():

    sampleList = []
    bCount = 0
    # If merging multiple bamfiles
    if args.i:
        print "Looking for files in directory: ", SAMdir
        print "Using the GTF file:", GTFFile
        print "These are defined in the script header. Can be changed permanently there, or for one offs using -D and -G"

        samlist = args.i.split(',')
        i = 1
        for file in samlist:
            print "Reading in file:", args.d + file
            if args.d:
                filepath = args.d + file
                f = pysam.AlignmentFile(filepath, "rb")
            else:
                f = pysam.AlignmentFile(file, "rb")
            sampleList.append(f)
            i += 1
    else:
        print "Need an input file! or use -h to see options"
        quit()

    # Get TSS/exon/UTR positions for gene model
    geneID = args.g

    if args.b:
        bedPosits = open(args.b, 'rb')
        for bedLine in bedPosits:
            bCount +=1

    if args.p == "2D":
        globalX = []
        globalY = []

        # For single region extract
        if args.b == False:
            # Define regions for extraction
            boundD = defineRegions()

            # Force min/max
            globalX.append(boundD['rangeStart'])
            globalX.append(boundD['rangeEnd'])
            globalY.append(args.s)
            globalY.append(args.S)

            for samfile in sampleList:
                print "Processing samfile"
                simpleExtract(globalX,globalY,samfile,boundD)
                samfile.close()

        # For multi region extract
        if args.b:
            # Force min/max
            globalX.append(-args.e)
            globalX.append(args.e)
            globalY.append(51)
            globalY.append(args.S)

            for samfile in sampleList:
                print "Processing samfile"
                # Doesn't need boundD as not opened bedfile yet
                multiRegionSimpleExtract(globalX,globalY,samfile)
                samfile.close()

        make2dHist(globalX,globalY,bCount)

    elif args.p == "3D":
        if len(samlist) >1:
            print "Currently not able to take multiple input files. Only use one at a time"

        boundD = defineRegions()

        for samfile in sampleList:
            print "Processing samfile"
            surfaceMat = regionExtract(samfile, boundD)
            samfile.close()

        make3dPlot(surfaceMat, boundD)

    print "Runtime:" + str(datetime.now() - startTime)

def defineRegions():

    # Set dictionary for all boundaries
    boundD = defaultdict(dict)

    # Load start/end regions
    if args.t is True:
        print "Using built in test co-ordinates: Chr1 16869721 16873927 | ~!NOT REAL RUN DATA!~ |"
        boundD['chrom'] = 'Chr1'
        boundD['TSS'] = 16869721
        boundD['TTS'] = 16873927

    if args.G is None:
        print "Using default gtf file location as none was defined:", args.G

    elif args.g is not None:
        gtfFile = open(args.G, 'rb')
        gene = args.g

        UTR5match = False
        UTR3match = False

        for line in gtfFile:
            if re.search(gene, line):
                if UTR5match == False and re.search("5UTR", line):
                    gline = line.split('\t')
                    boundD['chrom'] = gline[0]
                    if gline[6] == "+":
                        boundD['TSS'] = int(gline[3])
                    else:
                        boundD['TSS'] = int(gline[4])
                    UTR5match = True

                if UTR3match == False and re.search("3UTR", line):
                    gline = line.split('\t')
                    if gline[6] == "+":
                        boundD['TTS'] = int(gline[4])
                    else:
                        boundD['TTS'] = int(gline[3])
                    UTR3match = True

        print "Gene Boundaries:"
        print boundD['chrom'],min(boundD['TSS'],boundD['TTS']),max(boundD['TSS'],boundD['TTS'])

    elif args.b is None and args.B is None:
        print "Exiting as no co-ordinates given. Provide a gene using -g"
        quit()

    elif args.b is not None:
        print "Assessing multiple positions, in a 0 to ", args.e, "range"

    elif args.B is not None:
        print "Setting extraction co-ordinates from command line"

        b1 = args.B.split(':')
        boundD['chrom'] = b1[0]
        b2 = b1[1].split('-')
        boundD['TSS'] = int(b2[0])
        boundD['TTS'] = int(b2[1])

    boundD['glen'] = int(boundD['TSS'] - boundD['TTS'])
    boundD['rangeStart'] = min(boundD['TSS'],boundD['TTS']) - args.e
    boundD['rangeEnd'] = max(boundD['TSS'],boundD['TTS']) + args.e

    return boundD

def regionExtract(samfile,boundD):

    sizeDict = defaultdict(dict)

    x = 0
    while x <= 1001:
        y = boundD['rangeStart']
        while y <= boundD['rangeEnd']:
            sizeDict[x][y] = 0
            #print x, y
            y +=  args.r
        x += args.q

    sizeDict[0][boundD['rangeStart']] = 10000

    for read in samfile.fetch(boundD['chrom'],boundD['rangeStart'],boundD['rangeEnd']):
        if (abs(read.template_length) < args.s) and (abs(read.template_length) < args.S):

            pos = (int(read.reference_start / args.r)) * args.r
            size = (int(abs(read.template_length) / args.q)) * args.q

            while pos <= read.reference_end:
                if pos in sizeDict[size]:
                    sizeDict[size][pos] += 1
                #else:
                #    sizeDict[size][pos] = 1
                pos += 1


    # Sort and save dictionary by size and position
    surfaceMat = []

    for k in sorted(sizeDict.keys()):
        internalList = []
        for j in sorted(sizeDict[k].keys()):
            internalList.append(np.log(sizeDict[k][j] + 1))
        surfaceMat.append(list(internalList))

    return surfaceMat

def simpleExtract(globalX,globalY,samfile,boundD):
    leftpoint = min(boundD['rangeStart'],boundD['rangeEnd'])
    rightpoint = max(boundD['rangeStart'],boundD['rangeEnd'])

    for read in samfile.fetch(boundD['chrom'],leftpoint,rightpoint):
        if (abs(read.template_length) > args.s) and (abs(read.template_length) < args.S):

            pos = (int(read.reference_start / args.r)) * args.r
            size = (int(abs(read.template_length) / args.q)) * args.q

            while pos <= read.reference_end:
                globalX.append(pos)
                globalY.append(size)
                pos += 1

def multiRegionSimpleExtract(globalX,globalY,samfile):
    print "Reading bedfile"
    with open(args.b, 'rb') as f:
        for line in f:
            # take just first block before tab incase bed is followed by other annotation e.g. strand info
            bed0 = line.split('\t')
            # Chrom data
            b1 = bed0[0].split(':')
            chrom = b1[0]
            b2 = b1[1].split('-')
            midpoint = (int((int(b2[0]) + int(b2[1])) / 2) / args.r) * args.r
            regionRange = [midpoint - args.e, midpoint + args.e]

            for read in samfile.fetch(chrom,regionRange[0],regionRange[1]):
                if (abs(read.template_length) > args.s) and (abs(read.template_length) < args.S):

                    # Rounded matrix
                    pos = (int(abs(read.reference_start) / args.r)) * args.r
                    size = (int(abs(read.template_length) / args.q)) * args.q

                    while pos <= read.reference_end:
                        globalX.append(pos - midpoint)
                        globalY.append(size)
                        pos += args.r

def buildGeneModel():
    gtfFile = open(args.G, 'rb')
    gene = args.g
    geneModel = []
    gmin = None
    gmax = None

    for line in gtfFile:
        if re.search(gene, line):
            gline = line.split('\t')
            if gmin is None:
                gmin = gline[3]
                gmax = gline[3]

            # Find min
            if gline[3] < gmin:
                gmin = gline[3]
            elif gline[3] > gmax:
                gmax = gline[3]

            # Find max
            if gline[4] < gmin:
                gmin = gline[4]
            elif gline[4] > gmax:
                gmax = gline[4]

            if re.search('UTR', gline[2]):
                geneModel.append(
                    {
                    'type': 'rect',
                    'x0':gline[3] ,'y0': 11,'x1':gline[4] ,'y1': 19,
                    'line': {'color': 'rgb(0,0,0)','width': 2},
                    'fillcolor': 'rgb(0,100,0)'
                    }
                )
            if re.match('CDS', gline[2]):
                geneModel.append(
                    {
                    'type': 'rect',
                    'x0':gline[3] ,'y0': 11,'x1':gline[4] ,'y1': 19,
                    'line': {'color': 'rgb(0,0,0)','width': 2},
                    'fillcolor': 'rgb(0,0,0)'
                    }
                )
    # Make full transcript line
    geneModel.insert(0,
        {
        'type': 'rect',
        'x0':gmin ,'y0': 12.5,'x1':gmax ,'y1': 17.5,
        'line': {'color': 'rgb(0,0,0)','width': 2},
        'fillcolor': 'rgb(100,0,0)'
        }
    )

    return geneModel

def make3dPlot(surfaceMatrix, boundD):
    upstream = (float(args.e) / boundD['glen']) * 100

    # position in [] is x co-ord, position in [[]] is y co-ord, number is z co-ord
    data = [go.Surface(
        z=surfaceMatrix,
        zmin=0,
        zmax=5,
        colorscale=[[0, 'rgb(255,255,255)'], [0.01, 'rgb(50,50,50)'], [0.45, 'rgb(178,223,138)'], [1, 'rgb(227,26,28)']],
        contours=dict(),
        opacity=1
    )]

    layout = go.Layout(
        title=str(args.i),
        autosize=True,
        width=1590,
        height=1000,

        scene=dict(
            aspectratio=dict(x=4, y=2, z=1),
            aspectmode = 'manual',
            yaxis=dict(
                title='Particle Size',
                titlefont=dict(
                    family='Arial, sans-serif',
                    size=18,
                    #color='lightgrey'
                ),
                # log 50 and log 500. Christ knows why it wont let me log it in place...
                range=[1.698970004,2.698970004],
                type='log',
                autorange=True,
            ),
            xaxis=dict(
                title="Genomic Co-ordinates",
                titlefont=dict(
                    family='Arial, sans-serif',
                    size=31,
                    #color='lightgrey'
                ),
                #autorange=False,
                ticks='inside',
                showgrid=True,
                ticktext=["-" + str(args.e),str(boundD['rangeStart'] + args.e), str(boundD['rangeEnd'] - args.e), "+" + str(args.e)],
                tickvals=[0, upstream, (100 - upstream), 100]
            ),
            #zaxis=dict(range=[0,10])
        ),

    )
    #shapes = buildGeneModel()

    fig = go.Figure(data=data, layout=layout)
    outname = str(args.i) + "-3Dsurface.html"
    plotly.offline.plot(fig, filename=outname)

def make2dHist(poslist, sizelist,bCount):

    if args.g:
        geneModel = buildGeneModel()
        chartTitle = str("Gene:" + args.g + "  |||| Inputs:" + args.i)
        outname = str(args.g) + "-2Dcontour.html"
    elif args.B:
        geneModel = []
        chartTitle = str("Region: " + chrom + ":" + str(min(poslist)) + "-" + str(max(poslist)) + "|||| Inputs: " + args.i)
        outname = str(chrom + ":" + min(poslist) + "-" + max(poslist)) + "-2Dcontour.html"
    else:
        geneModel = []
        chartTitle = str("Inputs: " + args.i + " |||| No. Regions: " + str(bCount))
        outname = str(args.b) + "-2Dcontour.html"

    trace1 = go.Histogram2dContour(
        x=poslist,
        y=sizelist,
        zmin=0,
        zmax= args.z,
        contours=dict(showlines=False),
        colorscale=[[0, 'rgb(255,255,255)'], [0.25, 'rgb(31,120,180)'], [0.45, 'rgb(178,223,138)'], [0.65, 'rgb(51,159,44)'], [0.85, 'rgb(251,154,153)'], [1, 'rgb(227,26,28)']],
    )

    data = go.Data([trace1])

    layout = go.Layout(
        title=chartTitle,
        autosize=True,
        width=1800,
        height=800,
        yaxis=dict(
            title='Particle Size',
            titlefont=dict(
                family='Arial, sans-serif',
                size=18,
                #color='lightgrey'
            ),
            showgrid=True,
            type='log',
            #autorange=True
            # log 50 and log 500. Christ knows why it wont let me log it in place...
            range=[1,2.698970004]
        ),
        xaxis=dict(
            title="Genomic Co-ordinates",
            titlefont=dict(
                family='Arial, sans-serif',
                size=18,
                #color='lightgrey'
            )
        ),
        shapes = geneModel
    )

    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=outname)

if __name__ == "__main__":
    main()
