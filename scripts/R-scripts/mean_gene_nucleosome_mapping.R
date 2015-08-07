library(ggplot2)
library(reshape2)
library(scales)
library(vegan)

l <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-standard/profile_TSS_heatmap/genes/light-sigup.tsv", header=TRUE, check.names = FALSE)
d <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-standard/profile_TSS_heatmap/genes/dark-sigup.tsv", header=TRUE, check.names = FALSE)

l.annots <- l[1:4]
l[,1] <- NULL           ### DO 4 TIMES UNTIL I WORK OUT HOW TO REPEAT!

l.dec <- decostand(l, 'range', MARGIN=1)
l2 <- cbind(l.annots,l.dec)
l.melt <-melt(l2, id=c("Name","Exposure", "Treatment", "Gene"))
l.means <- dcast(l.melt, Name + Exposure + Treatment ~ variable, mean)
l.meanmelt <-melt(l.means, id=c("Name","Exposure", "Treatment"))

d.annots <- d[1:4]
d[,1] <- NULL            ### DO 4 TIMES UNTIL I WORK OUT HOW TO REPEAT!

d.dec <- decostand(d, 'range', MARGIN=1)
d2 <- cbind(d.annots,d.dec)
d.melt <-melt(d2, id=c("Name","Exposure", "Treatment", "Gene"))
d.means <- dcast(d.melt, Name + Exposure +Treatment ~ variable, mean)
d.meanmelt <-melt(d.means, id=c("Name","Exposure", "Treatment"))


pl <- ggplot(data=l.meanmelt, aes(x=as.numeric(as.character(variable)), y=value, colour=Exposure))
pd <- ggplot(data=d.meanmelt, aes(x=as.numeric(as.character(variable)), y=value, colour=Exposure))

p1 <- pl + 
  #coord_cartesian(ylim=c(0,150)) +#, xlim=c(-,270)) +
  stat_smooth(method="loess", span=0.1, se=TRUE) + 
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_vline(x=0, colour="black", lty=2) +
  #geom_vline(x=168, colour="black", lty=2) +
  scale_colour_brewer(palette="Set1") +
  facet_wrap(~ Treatment, ncol=1) +
  labs(title = "Genes under significantly higher expression during LIGHT exposure\n(Fischer exact test, p<0.05, FDR correction, >1CPM)")

p2 <- pd + 
  #coord_cartesian(ylim=c(0,150)) +#, xlim=c(-,270)) +
  stat_smooth(method="loess", span=0.1, se=TRUE) + 
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_vline(x=0, colour="black", lty=2) +
  #geom_vline(x=168, colour="black", lty=2) +
  scale_colour_brewer(palette="Set1") +
  facet_wrap(~ Treatment, ncol=1) +
  labs(title = "Genes under significantly higher expression during DARK exposure\n(Fischer exact test, p<0.05, FDR correction, >1CPM)")

multiplot(p1,p2, cols=2)
  


