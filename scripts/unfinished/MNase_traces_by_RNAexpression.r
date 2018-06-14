library(ggplot2)
library(data.table)
library(scales)
library(reshape2)
library(plyr)
library(vegan)


dataFiles <- lapply(Sys.glob("/home/sbi6dap/Projects/AGM/MNaseseq/particles/150_norm/WholeGenome/WholeGenome_TSS_heatmap/A*.xls"), read.csv)
RNAseq_data <- 
summary(dataFiles)
for (i in dataFiles){
  
}


tmp.melt <- melt(tmp.dec, id=c("pos"))






#######################################
p <-ggplot(data=tmp.melt, aes(x=as.numeric(as.character(pos)), y=value, colour=variable)) # BASIC
p +
  stat_smooth(method="loess", span=0.01, se=FALSE) +
  #geom_line(aes(x=as.numeric(as.character(variable)), y=value)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  scale_colour_brewer(palette="Paired") +
  geom_vline(xintercept=0, lty=2) +
  #scale_y_continuous(limits=c(-3.5, 3.5)) +
  theme(legend.position = "bottom") +
  labs(title = "TSS - lt120")