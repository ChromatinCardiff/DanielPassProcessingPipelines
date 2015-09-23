library(ggplot2)
library(data.table)
library(scales)
library(reshape2)
library(plyr)
library(vegan)

####################################################################################################################
###  BASIC CHART  ##################################################################################################
####################################################################################################################

# Read in
TSS <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/median/median_profile_TSS.xls", header=TRUE, sep="\t", row.names=1)
TTS <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/median/median_profile_TTS.xls", header=TRUE, sep="\t")
CSS <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/median/median_profile_CSS.xls", header=TRUE, sep="\t")
CTS <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/median/median_profile_CTS.xls", header=TRUE, sep="\t")
ESS <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/median/median_profile_ESS.xls", header=TRUE, sep="\t")
ETS <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/median/median_profile_ETS.xls", header=TRUE, sep="\t")
gene <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/median/median_profile_gene.xls", header=TRUE, sep="\t")

# Annotations
Exposure <- c("Light","Light","Dark","Dark","Light","Light","Dark","Dark")
Treatment <- c("Low","High","Low","High","Low","High","Low","High")

############### TSS ###############
# Processing
#rownames(TSS) <-TSS$pos
#TSS$pos <- NULL
TSS.dec <- as.data.frame(decostand(t(TSS), 'standardize', MARGIN=1))
TSS.dec$pos <- rownames(TSS.dec)
TSS.dec = cbind(Exposure,Treatment,TSS.dec)
#TSS.t <- t(TSS.dec)
TSS.melt <- melt(TSS.dec, id=c("pos","Exposure","Treatment"))

# Chart all columns
p <-ggplot(data=TSS.melt, aes(x=as.numeric(as.character(variable)), y=value, colour=Exposure))
TSS.plot <-p +
  stat_smooth(method="loess", span=0.01, se=TRUE) + 
  #geom_line(aes(x=as.numeric(as.character(variable)), y=value))+#, colour=Exposure)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  #scale_colour_brewer(palette="Paired") +
  geom_vline(x=0, colour="blue", lty=2) +
  #scale_y_continuous(limits=c(-3.5, 3.5)) +
  theme(legend.position = "right") +
  labs(title = "TSS") +
  facet_wrap(~ Treatment, ncol=1)
TSS.plot
############### TTS ###############
# Processing
rownames(TTS) <-TTS$pos
TTS.dec <- decostand(TTS, 'standardize', MARGIN=2)
TTS.dec$pos <- as.numeric(rownames(TTS))
TTS.melt <- melt(TTS.dec, id=c("pos"))

# Chart all columns
p <-ggplot(data=TTS.melt)
TTS.plot <-p +
  geom_line(aes(x=pos, y=value, colour=variable)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  scale_colour_brewer(palette="Paired") +
  geom_vline(x=0, colour="blue", lty=2) +
  #scale_y_continuous(limits=c(-3.5, 3.5)) +
  theme(legend.position = "none") +
  labs(title = "TTS")

############### CSS ###############
# Processing
rownames(CSS) <-CSS$pos
CSS.dec <- decostand(CSS, 'standardize', MARGIN=2)
CSS.dec$pos <- as.numeric(rownames(CSS))
CSS.melt <- melt(CSS.dec, id=c("pos"))

# Chart all columns
p <-ggplot(data=CSS.melt)
CSS.plot <-p +
  geom_line(aes(x=pos, y=value, colour=variable)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  scale_colour_brewer(palette="Paired") +
  geom_vline(x=0, colour="blue", lty=2) +
  #scale_y_continuous(limits=c(-3.5, 3.5)) +
  theme(legend.position = "none") +
  labs(title = "CSS")

############### CTS ###############
# Processing
rownames(CTS) <-CTS$pos
CTS.dec <- decostand(CTS, 'standardize', MARGIN=2)
CTS.dec$pos <- as.numeric(rownames(CTS))
CTS.melt <- melt(CTS.dec, id=c("pos"))

# Chart all columns
p <-ggplot(data=CTS.melt)
CTS.plot <-p +
  geom_line(aes(x=pos, y=value, colour=variable)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  scale_colour_brewer(palette="Paired") +
  geom_vline(x=0, colour="blue", lty=2) +
  #scale_y_continuous(limits=c(-3.5, 3.5)) +
  theme(legend.position = "none") +
  labs(title = "CTS")

############### ESS ###############
# Processing
rownames(ESS) <-ESS$pos
ESS.dec <- decostand(ESS, 'standardize', MARGIN=2)
ESS.dec$pos <- as.numeric(rownames(ESS))
ESS.melt <- melt(ESS.dec, id=c("pos"))

# Chart all columns
p <-ggplot(data=ESS.melt)
ESS.plot <-p +
  geom_line(aes(x=pos, y=value, colour=variable)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  scale_colour_brewer(palette="Paired") +
  geom_vline(x=0, colour="blue", lty=2) +
  #scale_y_continuous(limits=c(-3.5, 3.5)) +
  theme(legend.position = "none") +
  labs(title = "ESS")

############### ETS ###############
# Processing
rownames(ETS) <-ETS$pos
ETS.dec <- decostand(ETS, 'standardize', MARGIN=2)
ETS.dec$pos <- as.numeric(rownames(ETS))
ETS.melt <- melt(ETS.dec, id=c("pos"))

# Chart all columns
p <-ggplot(data=ETS.melt)
ETS.plot <-p +
  geom_line(aes(x=pos, y=value, colour=variable)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  scale_colour_brewer(palette="Paired") +
  geom_vline(x=0, colour="blue", lty=2) +
  #scale_y_continuous(limits=c(-3.5, 3.5)) +
  theme(legend.position = "none") +
  labs(title = "ETS")

############### gene ###############
# Processing
rownames(gene) <-gene$pos
gene.dec <- decostand(gene, 'standardize', MARGIN=2)
gene.dec$pos <- as.numeric(rownames(gene))
gene.melt <- melt(gene.dec, id=c("pos"))

# Chart all columns
p <-ggplot(data=gene.melt)
gene.plot <-p +
  geom_line(aes(x=pos, y=value, colour=variable)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  scale_colour_brewer(palette="Paired") +
  geom_vline(x=0, colour="blue", lty=2) +
  #scale_y_continuous(limits=c(-3.5, 3.5)) +
  theme(legend.position = "none") +
  labs(title = "gene")

multiplot(TSS.plot,TTS.plot,CSS.plot,CTS.plot,ESS.plot,ETS.plot, cols=3)
