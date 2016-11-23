library(ggplot2)
library(reshape)
library(scales)

# Genomic lengths for midlines
# Athal
chrom = c("Chr1","Chr2","Chr3","Chr4","Chr5")
chromlen <- c(30427680,19698290,23459840,18585060,26975510)
topArm = c(14449213,3607091,13590268,3052108,11132192)
bottomArm = c(14655898,16039854,9582349,14497759,14803217)
bottomMinusLen = chromlen - bottomArm

chromsizes <- data.frame(chrom, chromlen)
centros <- data.frame(chrom,topArm, bottomMinusLen)

# Data in
x <- read.table("/home/sbi6dap/Projects/NON-ATHAL/Yeast/prinegWT_tcv_1kb.sgr", sep="\t")
x <- read.table("/home/sbi6dap/Projects/NON-ATHAL/Yeast/nDNA_tcv_1kb.sgr", sep="\t")

colnames(x) <- c("chrom", "pos", "value")
#x$pos <- cbind(pos=row.names(x))

# Not needed if reading sgr
#x.melt <- melt(x, id="pos")
#summary(x.melt)

# Cut out <5% for better viewing
percentile5 <- quantile(x$value, .25, na.rm =TRUE)
x$transvalue <- ifelse((x$value-percentile5>0),x$value-percentile5,1)

prinegWT <- x
#x <- prinegWT
nDNA <- x

y <- data.frame(cbind(chrom=nDNA$chrom, pos=nDNA$pos, value=(prinegWT$value / nDNA$value)))

# Cut out <5% for better viewing
percentile <- quantile(y$value, .25, na.rm =TRUE)
y$transvalue <- ifelse((y$value-percentile>0),y$value-percentile,0)
head(y)

## Genomic plot
p <- ggplot(y, aes(x=as.numeric(pos),y=log10(transvalue)))
p + 
  theme_bw() +
  geom_ribbon(data=x, aes(ymin=(0 - log10(transvalue)), ymax = log10(transvalue))) +
  #coord_cartesian(xlim =  c(-2,2)) +
  scale_x_continuous(breaks = pretty_breaks(n=6)) +
  scale_y_continuous(breaks = pretty_breaks(n=3)) +
  #geom_segment(data=chromsizes, aes(x=0, y=0, xend=chromlen, yend=0), colour="blue") +
  #geom_segment(data=centros, aes(x=topArm, y=-1, xend=topArm, yend=1), colour="red", lty=1, size=1) +
  #geom_segment(data=centros, aes(x=bottomMinusLen, y=-1, xend=bottomMinusLen, yend=1), colour="red", lty=1, size=1) +
  facet_wrap(~chrom, nrow=1) + 
  coord_flip()
