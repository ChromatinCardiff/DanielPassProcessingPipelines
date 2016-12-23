library(ggplot2)
library(reshape)
library(scales)
library(vegan)

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
x <- read.table("/home/sbi6dap/Projects/AWT/wigs/tcv/ALD_tcv_10kb_merge.sgr", sep="\t")
x <- read.table("/home/sbi6dap/Projects/AWT/wigs/tcv/AtN_tcv_10kb.sgr", sep="\t")
#x <- read.table("/home/sbi6dap/Projects/NON-ATHAL/Yeast/prinegWT_tcv_1kb.sgr", sep="\t")
#x <- read.table("/home/sbi6dap/Projects/NON-ATHAL/Yeast/nDNA_tcv_1kb.sgr", sep="\t")

# gene density data # 
y <- read.table("/home/sbi6dap/Projects/REFDB/ARA_genedensity_10kb.txt", sep="\t")

head(x)
head(y)
xy <- cbind(x,y)

colnames(x) <- c("chrom", "pos", "value")

# Not needed if reading sgr
#x.melt <- melt(x, id="pos")
#summary(x.melt)

# Cut out <5% for better viewing
percentile5 <- quantile(x$value, .02, na.rm =TRUE)
x$transvalue <- ifelse((x$value-percentile5>0),x$value-percentile5,1)

# assign
ALD <- x
ATN <- x

# Not appropriate?
#ALD$stdvalue <- decostand(ALD$value, 'normalize', MARGIN=2)
#ATN$stdvalue <- decostand(ATN$value, 'normalize', MARGIN=2)

# Divide by number of reads * 1e9
ALD$stdvalue <- (ALD$value /  95806703.25) * 1000 * 1000 * 1000
ATN$stdvalue <- (ATN$value / 154394195) * 1000 * 1000 * 1000

# Divide by median
ALD$stdvalue <- ALD$value /  12.8
ATN$stdvalue <- ATN$value / 4.394

sub <- data.frame(chrom=ALD$chrom, pos=ALD$pos, 
                  ALDvalue=ALD$transvalue, 
                  ATNvalue=ATN$transvalue,
                  ALDstdvalue=ALD$stdvalue, 
                  ATNstdvalue=ATN$stdvalue, 
                  pcent=((ALD$stdvalue - ATN$stdvalue)/ALD$stdvalue)*100,
                  abspcent=ifelse(
                    ((ALD$stdvalue - ATN$stdvalue)/ALD$stdvalue)*100>0,
                    ((ALD$stdvalue - ATN$stdvalue)/ALD$stdvalue)*100,0)
                  )
head(sub)
summary(sub)

## Genomic plot
p <- ggplot(sub, aes(x=as.numeric(pos),y=log10(pcent)))
p + 
  theme_bw() +
  geom_bar(data=sub, stat="identity", aes(x=as.numeric(pos), y=5, width=150, colour=log10(abspcent+1), na.rm=TRUE)) +
  # Low = Naked occupancy > Experiment, High = Naked occupancy < Experiment #
  scale_colour_gradient(low = "white", high = muted("red")) +
  geom_ribbon(data=sub, aes(ymin=(0 - log10(ATNstdvalue)), ymax = log10(ATNstdvalue)), fill="green") +
  geom_ribbon(data=sub, aes(ymin=(0 - log10(ALDstdvalue)), ymax = log10(ALDstdvalue)), fill="black") +
  coord_cartesian(xlim =  c(-4,4)) +
  scale_x_continuous(breaks = pretty_breaks(n=6)) +
  scale_y_continuous(breaks = pretty_breaks(n=3)) +
  geom_segment(data=chromsizes, aes(x=0, y=0, xend=chromlen, yend=0), colour="blue") +
  geom_segment(data=centros, aes(x=topArm, y=-1, xend=topArm, yend=1), colour="red", lty=1, size=1) +
  geom_segment(data=centros, aes(x=bottomMinusLen, y=-1, xend=bottomMinusLen, yend=1), colour="red", lty=1, size=1) +
  facet_wrap(~chrom, nrow=1) +
  coord_flip()

multiplot(p1,p2, cols=1)
