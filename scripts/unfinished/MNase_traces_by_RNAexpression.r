library(ggplot2)
library(data.table)
library(scales)
library(reshape2)
library(plyr)
library(vegan)
library(ggrepel)

setwd("/home/sbi6dap/Projects/AGM/MNaseseq/particles/150_norm/TSS_for_R/")
filelist <- c("AGM07_150.Fnor.smooth.tsv","AGM08_150.Fnor.smooth.tsv","AGM09_150.Fnor.smooth.tsv","AGM10_150.Fnor.smooth.tsv","AGM11_150.Fnor.smooth.tsv","AGM12_150.Fnor.smooth.tsv")

AGM10 <- read.table("AGM10_150.Fnor.smooth.tsv", header=TRUE, sep="\t", row.names=1, check.names = FALSE)

RNAseq_data <- read.table("/home/sbi6dap/Projects/AGM/RNAseq/HTSeq/toptags_edgeR-WT_vs_cycd.csv", sep=",", header=TRUE)

summary(dataFiles)


files <- NULL
for (i in filelist){
  print(i)
  #fh <- substr(i,1,5)
  files[[substr(i,1,5)]] <- read.table(i, header=TRUE, sep="\t", row.names=1, check.names = FALSE)
}

a <- NULL
for (i in files){
  a[[deparse(substitute(i))]] <- cnc(i,6)
}


# cluster and chart function
cnc <- function(x,clusters){
  annot <- deparse(substitute(x))
  x.window <- x[120:160]
  k <- kmeans(x.window, clusters, iter.max=1000)
  x.window.k <- cbind(x, cluster=k$cluster)
  table(k$cluster)
  
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
    ggtitle(annot)
    #geom_label_repel(aes(label=clustertab))
  #facet_wrap(~cluster, ncol=1, scales = "free")
  #scale_y_continuous(limits=c(0, 200))
  
  return(x.clusters.plot)
  #y <- paste("~/Dropbox/Work/Projects/Athal-projects/AGM/kmeans/", deparse(substitute(x)), sep="")
  #ggsave(filename=y,x.clusters.plot)
}

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

