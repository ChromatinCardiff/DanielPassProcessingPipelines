#!/usr/bin/Rscript
# Source dpryan79
# https://github.com/dpryan79/Answers/blob/master/SEQanswers_42420/GTF2LengthGC.R
#source("http://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")
biocLite("rtracklayer")
biocLite("Rsamtools")

library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
setwd("~/Projects/REFDB/")

GTFfile = "Araport11_genes.20151202.gtf"
FASTAfile = "TAIR10_Chr.all.fasta"

#Load the annotation and reduce it
GTF <- import.gff(GTFfile, format="gtf", genome="UMD3.1.83", asRangedData=F, feature.type="exon")
grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementLengths(grl))

#Open the fasta file
FASTA <- FaFile(FASTAfile)
open(FASTA)

#Add the GC numbers
elementMetadata(reducedGTF)$nGCs <- letterFrequency(getSeq(FASTA, reducedGTF), "GC")[,1]
elementMetadata(reducedGTF)$widths <- width(reducedGTF)

#Create a list of the ensembl_id/GC/length
calc_GC_length <- function(x) {
    nGCs = sum(elementMetadata(x)$nGCs)
    width = sum(elementMetadata(x)$widths)
    c(width, nGCs/width)
}
LengthGC <-t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_GC_length))
colnames(LengthGC) <- c("Length", "GC")
LengthGC.df <- data.frame(LengthGC)
LengthGC.df$ID <- row.names(LengthGC.df)

write.table(LengthGC, file="GC_lengths.tsv", sep="\t")
