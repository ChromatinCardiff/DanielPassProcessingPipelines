library(ggplot2)
library(data.table)
library(scales)
library(reshape2)
library(plyr)
library(vegan)

#in
TSSA <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/mean_ARA11_GMC_TSS.xls", header=TRUE, sep="\t", row.names=1)
TSST <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/mean_TAIR/mean_profile_TSS.xls", header=TRUE, sep="\t", row.names=1)
CSSA <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/mean_ARA11/mean_CSS.xls", header=TRUE, sep="\t", row.names=1)
CSST <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/mean_TAIR/mean_profile_CSS.xls", header=TRUE, sep="\t", row.names=1)

# Annotations
Exposure <- c("Light","Light","Dark","Dark","Light","Light","Dark","Dark")
Treatment <- c("Low","High","Low","High","Low","High","Low","High")

############### TSS ###############

TSSA.dec <- as.data.frame(decostand(t(TSSA), 'standardize', MARGIN=1))
Annot <- c("ARA","ARA","ARA","ARA","ARA","ARA","ARA","ARA")
TSSA.dec$pos <- rownames(TSSA.dec)
TSSA.raw <- as.data.frame(t(TSSA))
TSSA.raw$pos <- rownames(TSSA.raw)
TSSA.raw = as.data.frame(cbind(Exposure,Treatment,Annot,TSSA.raw))
TSSA.dec = cbind(Exposure,Treatment,Annot,TSSA.dec)
TSSA.melt <- melt(TSSA.dec, id=c("pos","Exposure","Treatment", "Annot"))

TSST.dec <- as.data.frame(decostand(t(TSST), 'standardize', MARGIN=1))
Annot <- c("TAIR","TAIR","TAIR","TAIR","TAIR","TAIR","TAIR","TAIR")
TSST.raw <- as.data.frame(t(TSST))
TSST.raw$pos <- rownames(TSST.raw)
TSST.raw = as.data.frame(cbind(Exposure,Treatment,Annot,TSST.raw))
TSST.dec$pos <- rownames(TSST.dec)
TSST.dec = cbind(Exposure,Treatment,Annot,TSST.dec)
TSST.melt <- melt(TSST.dec, id=c("pos","Exposure","Treatment", "Annot"))

CSSA.dec <- as.data.frame(decostand(t(CSSA), 'standardize', MARGIN=1))
Annot <- c("ARA","ARA","ARA","ARA","ARA","ARA","ARA","ARA")
CSSA.raw = as.data.frame(cbind(Exposure,Treatment,Annot,t(CSSA)))
CSSA.raw$pos <- rownames(CSSA.raw)
CSSA.dec = cbind(Exposure,Treatment,Annot,CSSA.dec)
CSSA.melt <- melt(CSSA.raw, id=c("pos","Exposure","Treatment", "Annot"))

CSST.dec <- as.data.frame(decostand(t(CSST), 'standardize', MARGIN=1))
Annot <- c("TAIR","TAIR","TAIR","TAIR","TAIR","TAIR","TAIR","TAIR")
CSST.raw = as.data.frame(cbind(Exposure,Treatment,Annot,t(CSST)))
CSST.raw$pos <- rownames(CSST.raw)
CSST.dec = cbind(Exposure,Treatment,Annot,CSST.dec)
CSST.melt <- melt(CSST.raw, id=c("pos","Exposure","Treatment", "Annot"))

#combine melts
TSSall.melt <- rbind(TSST.melt,TSSA.melt)
CSSall.melt <- rbind(CSST.melt,CSSA.melt)

# Chart all columns
p <-ggplot(data=TSSall.melt, aes(x=as.numeric(as.character(variable)), y=value, colour=Treatment, lty=Annot))
TSS.plot <- p +
  stat_smooth(method="loess", span=0.1, se=TRUE) + 
  #geom_line(aes(x=as.numeric(as.character(variable)), y=value)) + #, colour=Exposure)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  #scale_colour_brewer(palette="Paired") +
  geom_vline(x=0, colour="blue", lty=2) +
  #scale_y_continuous(limits=c(-3.5, 3.5)) +
  theme(legend.position = "bottom") +
  labs(title = "TSS") +
  facet_wrap(~ Treatment, ncol=1)
TSS.plot

p <-ggplot(data=CSSall.melt, aes(x=as.numeric(as.character(variable)), y=value, colour=Exposure, lty=Annot))
CSS.plot <- p +
  stat_smooth(method="loess", span=0.01, se=FALSE) + 
  #geom_line(aes(x=as.numeric(as.character(variable)), y=value))+#, colour=Exposure)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  #scale_colour_brewer(palette="Paired") +
  geom_vline(x=0, colour="blue", lty=2) +
  #scale_y_continuous(limits=c(-3.5, 3.5)) +
  theme(legend.position = "bottom") +
  labs(title = "CSS") +
  facet_wrap(~ Treatment, ncol=1)
