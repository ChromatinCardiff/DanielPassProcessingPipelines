library(ggplot2)
library(reshape)
library(dplyr)
library(tibble)
#install.packages("ggsignif")
library(ggsignif)

IDtoGene <- read.delim("/home/daniel/GRH/YSR/Homo_sapiens.GRCh38.97.map.txt", header=T)



x <- data.frame(t(read.csv("/home/daniel/GRH/YSR/rmdup/temp.csv", sep="\t", row.names= 1)))
x2 <- x %>% rownames_to_column("Sample_ID")
meta <- read.csv("/home/daniel/GRH/YSR/targets.txt", sep="\t")
dat <- merge(x2, meta, by="Sample_ID")

dat.melt <-melt(dat)
dat.annot <- merge(dat.melt, IDtoGene, by.x="variable", by.y="Id" )

#t <- subset(dat.annot, Name=="HYDIN")


p <- ggplot(dat.annot, aes(Type, value))
p + geom_boxplot()  +
  facet_wrap(facets = "Name", scales="free")
