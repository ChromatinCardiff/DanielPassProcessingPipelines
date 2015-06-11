library(ggplot2)
library(reshape2)
library(scales)


i <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks/profile_TSS_heatmap/genes/tmp.txt", header=TRUE, check.names = FALSE)
row.names(i) <- i$row.names
t.i <- t(i)


i$Treatment <- row.names(i)

i.melt <-melt(i, id=c("Name","Exposure", "Treatment"))

p <- ggplot(data=i.melt, aes(x=as.numeric(as.character(variable)), y=value, colour=Treatment))
p + 
  #coord_cartesian(ylim=c(0,150)) +#, xlim=c(-,270)) +
  stat_smooth(method="loess", span=0.05, se=TRUE) + 
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_vline(x=0, colour="black", lty=2) +
  #geom_vline(x=168, colour="black", lty=2) +
  scale_colour_brewer(palette="Set1") +
  facet_wrap(~ Exposure, ncol=1)
  


