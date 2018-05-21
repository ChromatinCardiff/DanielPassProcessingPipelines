source("https://bioconductor.org/biocLite.R")
biocLite("karyoploteR")
library(karyoploteR)

a
kp <- plotKaryotype()

kp <- plotKaryotype(genome="hg19", plot.type=2, chromosomes=c("chr1", "chr2", "chr3"))