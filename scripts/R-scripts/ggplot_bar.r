library(ggplot2)
library(scales)
library(reshape2)

x <- read.table("Dropbox/Work/Manuscripts/Worm/Voxel/Data/L1raw.csv", header=TRUE, sep=",", check.names=FALSE)
x

x.melt <- melt(x, id=c("ID"))
head(x.melt)


p <- ggplot(x.melt, aes(x=ID, y = value))
p +
  geom_bar(stat = 'identity') +
  scale_x_continuous(breaks=seq(1,19))+
  coord_flip()+
  facet_wrap(~ variable, ncol=6)
