library(ggplot2)
library(data.table)
library(scales)
library(reshape2)
library(plyr)
library(vegan)

x <- read.table("/home/sbi6dap/Projects/ALD/totalcoverage/mean_ESS_heatmap/LDcommonmodel-ARA11.genepred.ess.ES09_totalcov.Fnor.wig.heatmap_1st-exon.xls.cut", header=TRUE, check.names=FALSE, sep="\t", row.names=1)
y <- read.table("/home/sbi6dap/Projects/ALD/totalcoverage/mean_ESS_heatmap/tmp.xls", header=TRUE, check.names=FALSE, sep="\t")
z <- read.table("/home/sbi6dap/Projects/ALD/totalcoverage/mean_ESS_heatmap/LDcommonmodel-ARA11.genepred.ess.ES09_totalcov.Fnor.wig.heatmap.xls.cut", header=TRUE, check.names=FALSE, sep="\t")

z$name <- NULL

x.mean <- colMeans(x)
y.mean <- colMeans(y)
z.mean <- colMeans(z)

x.dec <- decostand(x, 'standardize', MARGIN=2)
x.dec$pos <- as.numeric(rownames(x))


x.melt <- melt(x.mean)
x.melt$pos <- row.names(x.melt)
y.melt <- melt(y.mean)
y.melt$pos <- row.names(y.melt)
z.melt <- melt(z.mean)
z.melt$pos <- row.names(z.melt)

x.dec <- decostand(x.mean, method="range", 2)
y.dec <- decostand(y.mean, method="range", 2)
z.dec <- decostand(z.mean, method="range", 2)

x.decmelt <- melt(x.dec)
x.decmelt$pos <- row.names(x.melt)
y.decmelt <- melt(y.dec)
y.decmelt$pos <- row.names(y.melt)
z.decmelt <- melt(z.dec)
z.decmelt$pos <- row.names(z.melt)


p <- ggplot()
p1 <- p + 
  #geom_line(data=x.melt, aes(x=as.numeric(pos), y=value), col="red") +
  geom_line(data=y.melt, aes(x=as.numeric(pos), y=value), col="blue") +
  geom_line(data=z.melt, aes(x=as.numeric(pos), y=value), col="black") +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_vline(xintercept=0, colour="blue", lty=2) 
  
p2 <- p + 
  geom_line(data=x.melt, aes(x=as.numeric(pos), y=value), col="red") +
  #geom_line(data=y.melt, aes(x=as.numeric(pos), y=value), col="blue") +
  #geom_line(data=z.melt, aes(x=as.numeric(pos), y=value), col="black") +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_vline(xintercept=0, colour="blue", lty=2) +
  theme(legend.position = "right") 

multiplot(p1,p2)
