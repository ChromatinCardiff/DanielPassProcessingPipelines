#!/usr/bin/python

###################################################
## Daniel Pass | github.com/passdan | March 2016 ##
###################################################

# Basic processing
import argparse
import pysam
import re
from collections import defaultdict
# Numbers
import numpy as np
from scipy.stats.kde import gaussian_kde
from scipy import stats
from sklearn.neighbors import KernelDensity
# Graphing
import plotly
import plotly.graph_objs as go

## Define some parameters IF YOU WANT TO SAVE TYPING! ##
SAMdir = "/home/sbi6dap/RawData/AGM/MNase/bowtie/Col0-mapped-Lers/"
GTFFile = "/home/sbi6dap/Projects/REFDB/Araport11_genes.20151202.gtf"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Import Bam file and extract defined region for 2D or 3D plotting')
    parser.add_argument('-i', help='Input a Bam file or a list of comma-separated Bam files to merge (format: one.bam,two.bam,three.bam)')
    parser.add_argument('-d', default=SAMdir, help='Optional: Directory where bamfiles are (saves on typing!)')
    parser.add_argument('-o', help='Output prefix')
    parser.add_argument('-b', help='Regions to extract (Bed file)')
    parser.add_argument('-B', help='Individual set of co-ordinates in bed format (`Chr1:1000-1510`) (Requires surrounding quotes to escape semi-colon)')
    parser.add_argument('-g', help='Gene transcript to get positional information from | REQUIRES GTF FILE (-G)')
    parser.add_argument('-G', default=GTFFile, help='gtf file to get gene model information from | REQUIRES GENE/TRANSCRIPT ID')
    parser.add_argument('-e', default=500, type=int, help='Integer for upstream/downstream extension')
    parser.add_argument('-q', default=10, type=int, help='Integer for rounding counts (fragment size, default: 10)')
    parser.add_argument('-r', default=50, type=int, help='Integer for rounding counts (bases, default: 50)')
    parser.add_argument('-p', default="2D", help='Choose plot type (2D or 3D, default:2D)')
    parser.add_argument('-t', action='store_true', help='Use default co-ordinates for testing')

    args = parser.parse_args()

def main():

    sampleList = []
    # If merging multiple bamfiles
    if args.i:
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

    # Define regions for extraction
    boundD = defineRegions()

    #kde = kdeTransform(poslist,sizelist)

    if args.p == "2D":
        globalX = []
        globalY = []

        # Force min/max
        globalX.append(boundD['rangeStart'])
        globalX.append(boundD['rangeEnd'])
        globalY.append(51)
        globalY.append(1000)

        for samfile in sampleList:
            print "Processing samfile"
            simpleExtract(globalX,globalY,samfile,boundD)
            samfile.close()

        make2dHist(globalX,globalY,boundD)

    elif args.p == "3D":
        print "Currently broken for multiple input files. Only use one at a time"
        for samfile in sampleList:
            print "Processing samfile"
            surfaceMat = regionExtract(samfile, boundD)
            samfile.close()

        make3dPlot(surfaceMat, boundD)

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

        for line in gtfFile:
            if re.search(gene, line):
                if re.search("5UTR", line):
                    gline = line.split('\t')
                    boundD['chrom'] = gline[0]
                    if gline[6] == "+":
                        boundD['TSS'] = int(gline[3])
                    else:
                        boundD['TSS'] = int(gline[4])
                if re.search("3UTR", line):
                    gline = line.split('\t')
                    if gline[6] == "+":
                        boundD['TTS'] = int(gline[4])
                    else:
                        boundD['TTS'] = int(gline[3])

        print "Gene Boundaries:"
        print boundD['chrom'],min(boundD['TSS'],boundD['TTS']),max(boundD['TSS'],boundD['TTS'])

    elif args.b is None and args.B is None:
        print "Exiting as no co-ordinates given provide a gene using -g"
        quit()

    elif args.b is not None:
        print "retired method, use gene name and gtf instead"
        quit()

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
    #return(chrom,regionStart,regionEnd,regionLength)

def regionExtract(samfile,boundD):

    #sizeDict = {}
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
        #pos = read.reference_start
        #size = abs(read.template_length)

        #pos = int(round(read.reference_start, -1))
        #size = int(round(abs(read.template_length), -1))

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
        #pos = read.reference_start
        #size = abs(read.template_length)

        pos = (int(read.reference_start / args.r)) * args.r
        size = (int(abs(read.template_length) / args.q)) * args.q

        while pos <= read.reference_end:
            globalX.append(pos)
            globalY.append(size)
            pos += 1

    #return (m1,m2)

def kdeTransform(poslist,sizelist):
    #data=numpy.random.multivariate_normal(mu,surfaceMat,1000)
    #values = data.T
    #kde = stats.gaussian_kde(values)

    #values = np.vstack([poslist,sizelist])
    print type(poslist)
    print type(sizelist)
    samp = np.vstack([poslist,sizelist])
    kde = gaussian_kde(samp)
    #print values

    #kde = KernelDensity(kernel='gaussian', bandwidth=0.08804).fit(values)

    print kde
    return kde

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
                    'x0':gline[3] ,'y0': 51,'x1':gline[4] ,'y1': 59,
                    'line': {'color': 'rgb(0,0,0)','width': 2},
                    'fillcolor': 'rgb(0,100,0)'
                    }
                )
            if re.match('CDS', gline[2]):
                geneModel.append(
                    {
                    'type': 'rect',
                    'x0':gline[3] ,'y0': 51,'x1':gline[4] ,'y1': 59,
                    'line': {'color': 'rgb(0,0,0)','width': 2},
                    'fillcolor': 'rgb(0,0,0)'
                    }
                )
    # Make full transcript line
    geneModel.insert(0,
        {
        'type': 'rect',
        'x0':gmin ,'y0': 52.5,'x1':gmax ,'y1': 57.5,
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
        colorscale=[[0, 'rgb(255,255,255)'], [0.01, 'rgb(50,50,50)'], [0.45, 'rgb(178,223,138)'], [1, 'rgb(227,26,28)']],
        contours=dict(
            # x=dict(
            #     show=True,
            #     highlightwidth=1
            # ),
            # y=dict(
            #     show=True,
            #     highlightwidth=1
            # )
        ),
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
            zaxis=dict(range=[4,10])
        ),

    )
    #shapes = buildGeneModel()

    fig = go.Figure(data=data, layout=layout)
    outname = str(args.i) + "-3Dsurface.html"
    plotly.offline.plot(fig, filename=outname)

def make2dHist(poslist, sizelist, boundD):

    trace1 = go.Histogram2dContour(
        x=poslist,
        y=sizelist,
        zmin=0,
        zmax=10000,
        contours=dict(showlines=False),
        colorscale=[[0, 'rgb(255,255,255)'], [0.25, 'rgb(31,120,180)'], [0.45, 'rgb(178,223,138)'], [0.65, 'rgb(51,159,44)'], [0.85, 'rgb(251,154,153)'], [1, 'rgb(227,26,28)']],
    )

    data = go.Data([trace1])

    layout = go.Layout(
        title= str("Gene:" + args.g + "  |||| Inputs:" + args.i),
        #title= str("Inputs:" + args.i),
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
            range=[1.698970004,2.698970004]
        ),
        xaxis=dict(
            title="Genomic Co-ordinates",
            titlefont=dict(
                family='Arial, sans-serif',
                size=18,
                #color='lightgrey'
            )
        ),
        shapes = buildGeneModel()
    )

    fig = go.Figure(data=data, layout=layout)
    outname = str(args.o) + "-2Dcontour.html"
    plotly.offline.plot(fig, filename=outname)

if __name__ == "__main__":
    main()
