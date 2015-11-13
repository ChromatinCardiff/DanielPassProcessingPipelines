library(ggplot2)
library(vegan)
library(reshape2)
library(scales)
library(assertthat)
library(plyr)

x <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/mean_TAIR/mean_profile_TSS_heatmap/ES10-spline-avg.txt", header=TRUE, check.names=FALSE, sep="\t")
y <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/mean_TAIR/mean_profile_TSS_heatmap/all-1nucl-align.txt", header=TRUE, check.names=FALSE, sep="\t")


winsorize <-
  function(x, q=0.01)
  {
    assert_that(is.numeric(x))
    assert_that(is.number(q), q>=0, q<=1)
    
    lohi <- quantile(x, c(q, 1-q), na.rm=TRUE)
    if(diff(lohi) < 0) lohi <- rev(lohi)
    
    #x[!is.na(x) & x < lohi[1]] <- lohi[1]
    x[!is.na(x) & x > lohi[2]] <- lohi[2]
    x
  }

y.annots <- y[1:2]
x[,1] <-NULL
y[,1] <-NULL

y.win <- x
y.win[, -1] <- sapply(y.win[,-1], winsorize)
y.dec <- decostand(y.win, 'range', MARGIN=1)

## Averaging lines
y.dec <- cbind(y.annots, y.dec)
y.means <- aggregate(y.dec, list(Region = y.dec$sample), mean)
y.means$sample <- NULL
y.means$name <-NULL
y.melt <- melt(y.means, id=c("Sample", "Exposure", "Treatment"))
p <-ggplot(data=x.melt, aes(x=as.numeric(as.character(variable)), y=value))#, fill=Exposure))
p + 
  #stat_summary(geom="ribbon", fun.ymin="min", fun.ymax="max", alpha=0.5, colour="black", lty=5) +
  geom_line(aes(colour = Sample, lty=Exposure)) + 
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  #coord_cartesian(ylim = c(0.1,0.35), xlim = c(-125,125)) +  ##Standard
  coord_cartesian(ylim = c(0.1,0.35), xlim = c(-10,100))   ## +region
  #facet_wrap(~ Treatment, ncol=1)
  
  #scale_y_continuous(limits=c(0.1, 0.325)) 
  


write.table(x.means, file="/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/mean/mean_profile_TSS_heatmap/all-means-wind-dec-1nucl-align.txt")
write.table(y.means, file="/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/mean/mean_profile_TSS_heatmap/all-means-wind-dec-no-align.txt")
x.means <-read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/mean/mean_profile_TSS_heatmap/all-means-wind-dec-1nucl-align.txt", header=TRUE, check.names=FALSE)
y.means <-read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/mean/mean_profile_TSS_heatmap/all-means-wind-dec-no-align.txt", header=TRUE, check.names=FALSE)

#KMEANS
k <- kmeans(x.dec, 4)
x.k <- cbind(x.dec, cluster=k$cluster)
x.melt <- melt(x.k, id=c("cluster"))
x.means <- dcast(x.melt, cluster ~ variable, mean)
x.meanmelt <-melt(x.means, id=c("cluster"))

p <-ggplot(data=x.meanmelt, aes(x=as.numeric(as.character(variable)), y=value, colour=cluster))
p + geom_line()


#basic plot
#x <- x[1000:2781]
x.melt <- melt(x)
y.melt <- melt(y)
yES10.melt <- subset(y.melt, grepl("ES10", name))
p <-ggplot()
p + 
  geom_line(data=x.melt, aes(x=as.numeric(as.character(variable)), y=value),colour="red") +
  geom_line(data=y.melt, aes(x=as.numeric(as.character(variable)), y=value)) 
  stat_smooth(data=x.melt, aes(x=as.numeric(as.character(variable)), y=value), method="loess", span=0.01, se=FALSE) +
  stat_smooth(data=yES10.melt, aes(x=as.numeric(as.character(variable)), y=value, colour=name), method="loess", span=0.01, se=FALSE) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) 
  legend("topright")
  

