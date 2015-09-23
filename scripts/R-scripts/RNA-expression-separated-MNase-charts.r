library(ggplot2)
library(data.table)
library(assertthat)
library(scales)
library(reshape2)
library(plyr)
library(vegan)

####################################################################################################################
###  BASIC CHART  ##################################################################################################
####################################################################################################################

# Read in
x <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/diif-expression/light-primary.csv", header=TRUE, sep=" ")
summary(x)

# Processing
rownames(x) <-x$pos
x.win <- x
x.win[, -1] <- sapply(x.win[,-1], winsorize)
x.dec <- decostand(x.win, 'standardize', MARGIN=2)
x.dec$pos <- as.numeric(rownames(x))
x.melt <- melt(x.dec, id=c("pos"))

# Chart all columns
p2 <-ggplot(data=x.melt)
p2 +
  geom_line(aes(x=pos, y=value, colour=variable)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  scale_colour_brewer(palette="Paired") +
  geom_vline(x=0, colour="blue") 

# Separate x.melt into seperate dataframes
x.melt25 <- subset(x.melt, grepl("primary25", variable))
x.melt50 <- subset(x.melt, grepl("primary50", variable))
x.melt75 <- subset(x.melt, grepl("primary75", variable))
x.melt100 <- subset(x.melt, grepl("primary100", variable))

x.melt25_2k <- x.melt25[sample(1:nrow(x.melt25), 2000, replace=FALSE),]
x.melt50_2k <- x.melt50[sample(1:nrow(x.melt50), 2000, replace=FALSE),]
x.melt75_2k <- x.melt75[sample(1:nrow(x.melt75), 2000, replace=FALSE),]
x.melt100_2k <- x.melt100[sample(1:nrow(x.melt100), 2000, replace=FALSE),]

# Select columns
p1 <- ggplot(data=x.melt25) +
  geom_line(aes(x=pos, y=value, colour=variable)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_vline(x=0, colour="blue", lty=2) +
  labs(title = "Lowest 25% expressed genes (0% to 25%)", x = "Distance from TSS", y = "Normalised abundance (standardize)") +
  theme(legend.position="none")
p2 <- ggplot(data=x.melt50) +
  geom_line(aes(x=pos, y=value, colour=variable)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_vline(x=0, colour="blue", lty=2) +
  labs(title = "Medium-low expressed genes (25% to 50%)", x = "Distance from TSS", y = "Normalised abundance (standardize)") +
  theme(legend.position="none")
p3 <- ggplot(data=x.melt75) +
  geom_line(aes(x=pos, y=value, colour=variable)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_vline(x=0, colour="blue", lty=2) +
  labs(title = "Medium-high expressed genes (50% to 75%)", x = "Distance from TSS", y = "Normalised abundance (standardize)") +
  theme(legend.position="none")
p4 <- ggplot(data=x.melt100) +
  geom_line(aes(x=pos, y=value, colour=variable)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_vline(x=0, colour="blue", lty=2) +
  labs(title = "Highest 25% expressed genes (75% to 100%)", x = "Distance from TSS", y = "Normalised abundance (standardize)") +
  theme(legend.position="none")

multiplot(p1,p3,p2,p4, cols=2)


