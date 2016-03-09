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

#other
TSSA <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/mean_ARA11/mean_CSS.xls", header=TRUE, sep="\t", row.names=1)
TSST <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/mean_TAIR/mean_profile_CSS.xls", header=TRUE, sep="\t", row.names=1)

#tmp
tmp <- read.table("/home/sbi6dap/Projects/ALD/totalcoverage/mean_ESS_heatmap/LDcommonmodel-ARA11.genepred.ess.ES09_totalcov.Fnor.wig.heatmap_1st-exon.xls.cut", header=TRUE, sep="\t", row.names=1)

# Annotations
Exposure <- c("Light","Light","Dark","Dark","Light","Light","Dark","Dark")
Treatment <- c("Low","High","Low","High","Low","High","Low","High")

############### TSS ###############

TSSA.dec <- as.data.frame(decostand(t(TSSA), 'standardize', MARGIN=1))
Annot <- c("ARA","ARA","ARA","ARA","ARA","ARA","ARA","ARA")
TSSA.dec$pos <- rownames(TSSA.dec)
TSSA.dec = cbind(Exposure,Treatment,Annot,TSSA.dec)
TSSA.melt <- melt(TSSA.dec, id=c("pos","Exposure","Treatment", "Annot"))

TSST.dec <- as.data.frame(decostand(t(TSST), 'standardize', MARGIN=1))
Annot <- c("TAIR","TAIR","TAIR","TAIR","TAIR","TAIR","TAIR","TAIR")
TSST.dec$pos <- rownames(TSST.dec)
TSST.dec = cbind(Exposure,Treatment,Annot,TSST.dec)
TSST.melt <- melt(TSST.dec, id=c("pos","Exposure","Treatment", "Annot"))

#combine melts
TSSall.melt <- rbind(TSST.melt,TSSA.melt)
ESSall.melt <- rbind(ESST.melt,ESSA.melt)

# Chart all columns
p <-ggplot(data=TSSA.melt, aes(x=as.numeric(as.character(variable)), y=value, colour=Exposure, lty=Annot))
TSS.plot <-p +
  stat_smooth(method="loess", span=0.01, se=FALSE) + 
  #geom_line(aes(x=as.numeric(as.character(variable)), y=value))+#, colour=Exposure)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  #scale_colour_brewer(palette="Paired") +
  geom_vline(x=0, colour="blue", lty=2) +
  #scale_y_continuous(limits=c(-3.5, 3.5)) +
  theme(legend.position = "bottom") +
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

########### One Off charts ###########
tmp.dec <- as.data.frame(decostand(tmp, 'standardize', MARGIN=2, na.rm =TRUE))
tmp.dec$pos <- rownames(tmp.dec)
#tmp.dec = cbind(Exposure,Treatment,tmp.dec)
tmp.melt <- melt(tmp.dec, id=c("pos"))
tmp.melt <- melt(tmp.dec, id=c("pos", "Exposure","Treatment"))

tmp$pos <- rownames(tmp)
tmp.melt <- melt(tmp, id=c("pos"))
#######################################
p <-ggplot(data=tmp.melt, aes(x=as.numeric(as.character(pos)), y=value, colour=variable)) # BASIC
#p <-ggplot(data=tmp.melt, aes(x=as.numeric(as.character(variable)), y=value, colour=Exposure, lty=Treatment))
p +
  stat_smooth(method="loess", span=0.01, se=FALSE) +
  #geom_line(aes(x=as.numeric(as.character(variable)), y=value)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  scale_colour_brewer(palette="Paired") +
  #geom_vline(0, colour="blue", lty=2) +
  #scale_y_continuous(limits=c(-3.5, 3.5)) +
  theme(legend.position = "bottom") +
  labs(title = "TSS")
