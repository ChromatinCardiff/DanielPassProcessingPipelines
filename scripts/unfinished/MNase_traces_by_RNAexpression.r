library(ggplot2)
library(data.table)
library(scales)
library(reshape2)
library(dplyr)
library(vegan)
library(ggrepel)
library(gridExtra)
library(pandas)
import pandas as pd



setwd("/home/sbi6dap/Projects/AGM/MNaseseq/particles/150_norm/TSS_for_R/")
filelist <- c("AGM07_150.Fnor.smooth.tsv","AGM08_150.Fnor.smooth.tsv","AGM09_150.Fnor.smooth.tsv","AGM10_150.Fnor.smooth.tsv","AGM11_150.Fnor.smooth.tsv","AGM12_150.Fnor.smooth.tsv")

#AGM10 <- read.table("AGM10_150.Fnor.smooth.tsv", header=TRUE, sep="\t", row.names=1, check.names = FALSE)

RNAseq_data <- read.table("/home/sbi6dap/Projects/AGM/RNAseq/HTSeq/toptags_edgeR-WT_vs_cycd.csv", sep=",", header=TRUE)
RNAseq_counts <- read.table("/home/sbi6dap/Projects/AGM/RNAseq/HTSeq/AGM-WTG54.count", sep="\t", header=TRUE)
RNAseq_counts$WTavg <- rowMeans(RNAseq_counts[2:4])
RNAseq_counts$G54avg <- rowMeans(RNAseq_counts[5:7])

summary(dataFiles)


files <- NULL
for (i in filelist){
  print(i)
  #fh <- substr(i,1,5)
  files[[substr(i,1,5)]] <- read.table(i, header=TRUE, sep="\t", row.names=1, check.names = FALSE)
}

files[[1]]$ID <- row.names(files[[1]])
files[[2]]$ID <- row.names(files[[2]])
files[[3]]$ID <- row.names(files[[3]])
files[[4]]$ID <- row.names(files[[4]])
files[[5]]$ID <- row.names(files[[5]])
files[[6]]$ID <- row.names(files[[6]])


WT.df <- bind_rows(list(files[[1]], files[[2]], files[[3]])) %>% group_by(ID) %>% summarise_all(funs(mean))
G54.df <- bind_rows(list(files[[4]], files[[5]], files[[6]])) %>% group_by(ID) %>% summarise_all(funs(mean))

setwd("/home/sbi6dap/Dropbox/Work/Projects/Athal-projects/AGM/kmeans/")
ggsave(filename="AGM-WT-RNA.pdf",cnc(WT.df))
#cnc(files[[1]])

sapply(files,cnc)

# cluster and chart function
cnc <- function(x){
  clusters <- 6
  annot <- deparse(substitute(x))
  print(annot)
  x.window <- x[120:160]
  k <- kmeans(x.window, clusters, iter.max=1000)
  x.window.k <- cbind(x, cluster=k$cluster)
  
  # Make cluster table and RNA counts
  kTable <- data.frame(k$cluster)
  kTable$ID <- row.names(kTable)
  x.RNA <- merge(RNAseq_data, kTable, by="ID")
  head(x.RNA)
  #print(aggregate(x.RNA$logFC, list(x.RNA$k.cluster), mean))
  #clustTable <- cbind(RNA_FC=aggregate(x.RNA$logFC, list(x.RNA$k.cluster), mean), clusters=table(k$cluster))
  #tmp <- cbind(RNA_abun=data.frame(aggregate(x.RNA$WTavg, list(x.RNA$k.cluster), mean)), clusters=table(k$cluster))
  #clustTable <- tmp[, c(1,2,4)]
  
  #cluster size only
  #clustTable <- data.frame(table(k$cluster))
  
  # aggreate groups
  x.agg <- aggregate(x.window.k, list(x.window.k$cluster), mean)
  x.agg.melt <- melt(x.agg[2:ncol(x.agg-1)], id=c("cluster"))
  
  # Plot graph
  p <-ggplot(data=x.agg.melt, aes(x=as.numeric(as.character(variable)), y=value, colour=as.factor(cluster))) # BASIC
  x.clusters.plot <- p +
    geom_smooth(method="loess", span=0.01, se=FALSE)+
    scale_x_continuous(breaks = pretty_breaks(n=12)) +
    scale_colour_brewer(palette="Paired") +
    geom_vline(xintercept=0, lty=2) +
    ggtitle(annot) +
    annotation_custom(tableGrob(clustTable),  xmin=60, xmax=150, ymin=100, ymax=120)

  return(x.clusters.plot)
  #y <- paste("~/Dropbox/Work/Projects/Athal-projects/AGM/kmeans/", deparse(substitute(x)), sep="")
  #ggsave(filename=y,x.clusters.plot)
}

cnc(files[[1]])

AGM01.tmp.ksub <- ddply(AGM01.tmp.k, .(cluster), subset, sample(seq_along(cluster)<=200))
AGM01.tmp.ksub2 <- na.omit(AGM01.tmp.ksub)

#summary(AGM01.tmp.ksub)
AGM01.tmp.sort <- AGM01.tmp.ksub2[order(AGM01.tmp.ksub2[,ncol(AGM01.tmp.ksub2)]),]
AGM01.tmp.sort <- cbind(AGM01.tmp.sort, "idsort"=1:nrow(AGM01.tmp.sort))
#head(AGM01.tmp.sort, 100)


## Output gene cluster names
AGM01.tmp.genes <- subset(AGM01.tmp.sort, , c(42:43))        ######### Loot at column number
rownames(AGM01.tmp.sort) <- AGM01.tmp.sort[,104]
AGM01.tmp.annots <- AGM01.tmp.k[103]
write.table(AGM01.tmp.annots, file = "/home/sbi6dap/Projects/AGM/MNaseseq/AGM07-lt120clusters.tAGM01.tmpt")

AGM01.tmp.melt <- melt(AGM01.tmp.sort, id=c("cluster","idsort"))
#summary(AGM01.tmp.melt)
head(AGM01.tmp.melt, 10)


AGM07.window <- WT.df[120:160]
k <- kmeans(AGM07.window, clusters, iter.max=1000)
kTable <- data.frame(k$cluster)
kTable$ID <- row.names(kTable)
AGM07.RNA <- merge(RNAseq_data, kTable, by="ID")
aggregate(AGM07.RNA$AGM07, list(AGM07.RNA$k.cluster), mean)

table(AGM07.RNA$logCPM)


