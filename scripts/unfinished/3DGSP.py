#!/usr/bin/python

###################################################
## Daniel Pass | github.com/passdan | March 2016 ##
###################################################

# Processing
import argparse
import pysam
import numpy as np
from scipy.stats.kde import gaussian_kde
from scipy import stats
from collections import defaultdict
from sklearn.neighbors import KernelDensity
# Graphing
import plotly
import plotly.graph_objs as go



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Import Bam file and extract defined region for 2D or 3D plotting')
    parser.add_argument('-i', required=True, help='Input Bam file')
    parser.add_argument('-o', help='Output prefix')
    parser.add_argument('-b', help='Regions to extract (Bed file)')
    parser.add_argument('-B', help='Individual set of co-ordinates in Chr1,1000,1500 format (include commas)')
    parser.add_argument('-e', default=0, type=int, help='Integer for upstream/downstream extension')
    parser.add_argument('-r', default=10, type=int, help='Integer for rounding counts (size and base position)')
    parser.add_argument('-p', default="3D", help='Choose plot type (2D or 3D, default:3D)')
    parser.add_argument('-t', action='store_true', help='Use default co-ordinates for testing')

    args = parser.parse_args()

def main():
    # Load sam/bamfile
    samfile = pysam.AlignmentFile(args.i, "rb")

    # Load start/end regions
    if args.t is True:
        print "Using built in test co-ordinates: Chr1 16869721 16873927 | ~!NOT REAL RUN DATA!~ |"
        chrom = 'Chr1'
        regionStart = 16869721 - args.e
        regionEnd = 16873927 + args.e

    elif args.b is None and args.B is None:
        print "Exiting as no co-ordinates given"
        quit()

    elif args.b is not None:
        print "Not supported yet"
        quit()

    elif args.B is not None:
        print "Setting extraction co-ordinates from command line"

        bedParams = args.B.split(',')
        chrom = bedParams[0]
        regionStart = int(bedParams[1]) - args.e
        regionEnd = int(bedParams[2]) + args.e

    #cdsStart = 16869921

    #kde = kdeTransform(poslist,sizelist)

    if args.p == "2D":
        (poslist,sizelist) = simpleExtract(samfile,chrom, regionStart, regionEnd)
        samfile.close()

        make2dHist(poslist,sizelist, regionStart, regionEnd)

    elif args.p == "3D":
        surfaceMat = regionExtract(samfile,chrom, regionStart, regionEnd)
        samfile.close()

        make3dPlot(surfaceMat, regionStart, regionEnd)


def regionExtract(samfile,chr,start,end):

    #sizeDict = {}
    sizeDict = defaultdict(dict)

    x = 0
    while x <= 1001:
        y = start
        while y <= end:
            sizeDict[x][y] = 0
            #print x, y
            y +=  args.r
        x += args.r

    for read in samfile.fetch(chr, start, end):
        #pos = read.reference_start
        #size = abs(read.template_length)

        #pos = int(round(read.reference_start, -1))
        #size = int(round(abs(read.template_length), -1))

        pos = (int(read.reference_start / args.r)) * args.r
        size = (int(abs(read.template_length) / args.r)) * args.r

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
            internalList.append(sizeDict[k][j])
        surfaceMat.append(list(internalList))

    return surfaceMat

def simpleExtract(samfile,chr,start,end):
    m1 = []
    m2 = []

    # Force min/max
    m1.append(start)
    m1.append(start)
    m2.append(0)
    m2.append(1000)

    for read in samfile.fetch(chr, start, end):
        pos = read.reference_start
        size = abs(read.template_length)

        #pos = int(round(read.reference_start, -1))
        #size = int(round(abs(read.template_length), -1))

        #pos = (int(read.reference_start / args.r)) * args.r
        #size = (int(abs(read.template_length) / args.r)) * args.r

        while pos <= read.reference_end:
            m1.append(pos)
            m2.append(size)
            pos += 1

    return (m1,m2)

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


def make3dPlot(surfaceMatrix, regionStart, regionEnd):

    # position in [] is x co-ord, position in [[]] is y co-ord, number is z co-ord
    data = [go.Surface(z=surfaceMatrix)]
    layout = go.Layout(
        title=args.i,
        autosize=True,
        width=1600,
        height=1000,

        scene=dict(
            aspectratio=dict(x=3, y=1, z=1),
            aspectmode = 'manual',
            yaxis=dict(
                title='Particle Size',
                titlefont=dict(
                    family='Arial, sans-serif',
                    size=18,
                    #color='lightgrey'
                ),
                type='log',
                autorange=True,
            ),
            xaxis=dict(
                title="Genomic Co-ordinates",
                titlefont=dict(
                    family='Arial, sans-serif',
                    size=18,
                    #color='lightgrey'
                ),
                #type='log',
                #autorange=False,
                ticktext=[str(regionStart), str(regionEnd)],
                tickvals=[0,400]
            )
        )
    )

    fig = go.Figure(data=data, layout=layout)
    outname = str(args.i) + "-3Dsurface.html"
    plotly.offline.plot(fig, filename=outname)

def make2dHist(poslist, sizelist, regionStart, regionEnd):

    trace1 = go.Histogram2dContour(
        x=poslist,
        y=sizelist,
        contours=dict(
            showlines=False
        )
    )
    data= go.Data([trace1])

    layout = go.Layout(
        title=args.i,
        autosize=True,
        width=1600,
        height=600,
        yaxis=dict(
            title='Particle Size',
            titlefont=dict(
                family='Arial, sans-serif',
                size=18,
                #color='lightgrey'
            ),
            type='log',
            range=[50,1000]
        ),
        xaxis=dict(
            title="Genomic Co-ordinates",
            titlefont=dict(
                family='Arial, sans-serif',
                size=18,
                #color='lightgrey'
            )
        ),
        shapes = [{
            'type': 'line',
            'x0': 16869721,
            'y0': 50,
            'x1': 16873927,
            'y1': 50,
            'line': {
                'color': 'rgb(50, 171, 96)',
                'width': 10,
                'dash': 'dashdot',
            }
        }]
    )

    fig = go.Figure(data=data, layout=layout)
    outname = str(args.i) + "-2Dcontour.html"
    plotly.offline.plot(fig, filename=outname)

if __name__ == "__main__":
    main()
