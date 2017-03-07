library(ggplot2)
library(scales)
library(reshape2)
library(plyr)
library(vegan)
library(format)

####################################################################################################################
###  small_particles ###############################################################################################
####################################################################################################################

# Read in
setwd(/home/sbi6dap/Projects/AGM/MNaseseq/Col0/E2FC_targets/lt120traces/)
x.head <- read.table("/home/sbi6dap/Projects/AGM/MNaseseq/Col0/E2FC_targets/HEADER.txt", header=FALSE, sep="\t")

### Single file processing ###
x <- read.table("/home/sbi6dap/Projects/AGM/MNaseseq/Col0/E2FC_targets/lt120traces/trace-E2FC_target_list-L.xls.cut", header=FALSE, sep="\t", row.names=1)
colnames(x) <- x.head

WT <- x[grep("AGM07|AGM08|AGM09", rownames(x)),]
G54 <- x[grep("AGM10|AGM11|AGM12", rownames(x)),]

traces <- NULL
traces$WT <- apply(WT, 2, mean)
traces$G54 <- apply(G54, 2, mean)
traces.df <- data.frame(traces)
head(traces.df)
traces.df$pos <- row.names(traces.df)
traces.melt <- melt(traces.df, id=c("pos"))

L <- ggplot(data=traces.melt) +
  stat_smooth(aes(x=as.numeric(pos), y=value, colour=variable), method="loess", span=0.1, se=FALSE) +
  #geom_line(aes(x=as.numeric(pos), y=value, colour=variable)) +
  scale_x_continuous(breaks = pretty_breaks(n=6)) +
  geom_vline(xintercept=0, colour="blue", linetype="longdash") +
  labs(x = "Distance from TSS", y = "Abundance") +
  theme(legend.position="none")

}

### Multi File processing ###
letterstring <- "ABCDEFGHIJKL"
splitstring <- strsplit(letterstring,"")[[1]]

for (itter in splitstring){
  tmpfile <- paste("/home/sbi6dap/Projects/AGM/MNaseseq/Col0/E2FC_targets/Col0/lt120traces/trace-E2FC_target-", itter, ".xls.cut", sep="")
  x <- read.table(tmpfile, header=FALSE, sep="\t", row.names=1)
  colnames(x) <- x.head
  
  #WT <- x[grep("AGM07|AGM08|AGM09", rownames(x)),]
  #MOD <- x[grep("AGM10|AGM11|AGM12", rownames(x)),]
  
  WT <- x[grep("AGM01|AGM03|AGM05", rownames(x)),]
  MOD <- x[grep("AGM13|AGM14|AGM15", rownames(x)),]
  
  traces <- NULL
  traces$WT <- apply(WT, 2, mean)
  traces$MOD <- apply(MOD, 2, mean)
  traces.df <- data.frame(traces)
  traces.df$pos <- row.names(traces.df)
  head(traces.df)
  #traces.dec <- decostand(traces.df, 'standardize', MARGIN=2)
  #head(traces.dec)
  #traces.dec$pos <- row.names(traces.dec)
  traces.melt <- melt(traces.df, id=c("pos"))
  
  chartname <- toString(itter)
  assign(chartname, ggplot(data=traces.melt) +
    stat_smooth(aes(x=as.numeric(pos), y=value, colour=variable), method="loess", span=0.1, se=FALSE) +
    #geom_line(aes(x=as.numeric(pos), y=value, colour=variable)) +
    scale_x_continuous(breaks = pretty_breaks(n=6)) +
    geom_vline(xintercept=0, colour="blue", linetype="longdash") +
    labs(x = "TSS", y = "Abundance") +
    theme(legend.position="none")
  )
}


multiplot(B,C,D,E,F,G,H,I,J,K,L)
