library(ggplot2)
library(reshape)
library(miscTools)

sizes <- read.csv("~/Projects/ACS/analysis/dpos/profile_TSS_heatmap/gene_sizes.txt", header=TRUE)

ggplot(sizes, aes(len)) + 
  scale_y_reverse() + 
  stat_ecdf(aes(ymin=1, ymax=..y..),geom="ribbon") + 
  coord_cartesian(xlim = c(0, 5000))


temp <- read.csv("~/Downloads/temp.txt", header=FALSE)

x1 <- read.table("~/Projects/ALD/MNase-seq/dpos_peaks/profile_TSS_heatmap/uniq.ES09.txt", header=TRUE)
rownames(y1) <-x1$row.names
x1[,1] <- NULL
x2<-colMeans(x1, na.rm=TRUE)
x.melt <-melt(x2)
x.melt2 <-na.omit(x.melt)

y1 <- read.table("~/Projects/ALD/MNase-seq/dpos_peaks/profile_TSS_heatmap/masked.ES09.txt", header=TRUE)
rownames(y1) <-y1$row.names
y1[,1] <- NULL
y2<-colMeans(y1, na.rm=TRUE)
y.melt <-melt(y2)
y.melt2 <-na.omit(y.melt)

plot(y.melt$value, type="l", col="red", lty=1,  xaxt="n", xlab="", ylim=c(60,90), xlim=c(0,650))

plot(x.melt$value, type="l", col="black", lty=5,  xaxt="n", xlab="", ylim=c(60,90), xlim=c(0,650))
lines(y.melt$value, type="l", col="red")
abline(v=150, lty=1)
text(165, 60, "TSS")
abline(v=0, lty=1)
text(25, 60, "-1500bp")
abline(v=650, lty=1)
text(630, 60, "5000bp")
legend("topleft", c("Full Data", "Exons only"), fill=c("black", "red"))

mask <- read.csv("~/Projects/ALD/MNase-seq/dpos_peaks/profile_TSS_heatmap/genemask.csv", header=FALSE)
TSS.density <- colMode(mask)
plot(TSS.density, type="l")

apply(mask, 2, function(x) names(which.max(table(x))))

### dpos_profile map
profile <- read.table("~/Projects/ALD/MNase-seq/dpos_peaks/profile_TSS.xls", header=TRUE)
profile.melt <- melt(profile, id="pos")
p <- ggplot(data=profile.melt, aes(x=pos, y=value, colour=variable))
p + geom_line()


