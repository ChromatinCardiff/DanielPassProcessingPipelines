library(ggplot2)
library(reshape2)
library(scales)


i <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-standard/profile_TSS_heatmap/genes/LD-comp-test.tsv", header=TRUE, check.names = FALSE)
row.names(i) <- i$row.names
t.i <- t(i)

i.melt <-melt(i, id=c("Name","Exposure", "Treatment", "Gene"))

p <- ggplot(data=i.melt, aes(x=as.numeric(as.character(variable)), y=value, colour=Exposure))
p + 
  #coord_cartesian(ylim=c(0,150)) +#, xlim=c(-,270)) +
  stat_smooth(method="loess", span=0.05, se=TRUE) + 
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_vline(x=0, colour="black", lty=2) +
  #geom_vline(x=168, colour="black", lty=2) +
  scale_colour_brewer(palette="Set1") +
  facet_wrap(~ Gene, ncol=1)
  

ES09 <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq-80/profile_TSS_heatmap/ES09_80-1nuclalign-avg.txt", header=TRUE, check.names = FALSE, sep="\t")
ES10 <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq-80/profile_TSS_heatmap/ES10_80-1nuclalign-avg.txt", header=TRUE, check.names = FALSE, sep="\t")
ES11 <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq-80/profile_TSS_heatmap/ES11_80-1nuclalign-avg.txt", header=TRUE, check.names = FALSE, sep="\t")
ES12 <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq-80/profile_TSS_heatmap/ES12_80-1nuclalign-avg.txt", header=TRUE, check.names = FALSE, sep="\t")
ES13 <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq-80/profile_TSS_heatmap/ES13_80-1nuclalign-avg.txt", header=TRUE, check.names = FALSE, sep="\t")
ES14 <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq-80/profile_TSS_heatmap/ES14_80-1nuclalign-avg.txt", header=TRUE, check.names = FALSE, sep="\t")
ES15 <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq-80/profile_TSS_heatmap/ES15_80-1nuclalign-avg.txt", header=TRUE, check.names = FALSE, sep="\t")
ES16 <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq-80/profile_TSS_heatmap/ES16_80-1nuclalign-avg.txt", header=TRUE, check.names = FALSE, sep="\t")
a <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq-80/profile_TSS_heatmap/tmp.txt", header=TRUE, check.names = FALSE, sep="\t")


a.melt <-melt(a, id=c("name"))

ES09.melt <-melt(ES09, id=c("name"))
ES10.melt <-melt(ES10, id=c("name"))
ES11.melt <-melt(ES11, id=c("name"))
ES12.melt <-melt(ES12, id=c("name"))
ES13.melt <-melt(ES13, id=c("name"))
ES14.melt <-melt(ES14, id=c("name"))
ES15.melt <-melt(ES15, id=c("name"))
ES16.melt <-melt(ES16, id=c("name"))


p <-ggplot()

p + 
  coord_cartesian(xlim=c(-150,0)) +
  geom_line(data=a.melt, aes(x=as.numeric(as.character(variable)), y=(value), colour=name)) 
  geom_line(data=ES10.melt, aes(x=as.numeric(as.character(variable)), y=(value), colour="blue")) +
  geom_line(data=ES11.melt, aes(x=as.numeric(as.character(variable)), y=(value), colour="green")) 
  #geom_line(data=i.melt, aes(x=as.numeric(as.character(variable)), y=value, colour="red"))
  #stat_smooth(method="loess", span=0.1, se=FALSE) + 
  #scale_x_continuous(breaks = pretty_breaks(n=12)) +
  #geom_vline(x=0, colour="black", lty=2) 
  
