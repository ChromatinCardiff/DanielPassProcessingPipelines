library(vegan)
library(matrixStats)

setwd("/home/sbi6dap/Projects/ALD/totalcoverage/size_comparisons/ALL_CSS")
# Chromatin per-gene profile files
myFiles <- list.files(pattern="*_lt120.wig.heatmap.xls.uniq")
filenames = c("ES09","ES10","ES11","ES12","ES13","ES14","ES15","ES16")

for(x in seq_along(myFiles)){
 data.in <- read.table(myFiles[x], header=TRUE, sep="\t", row.names=1, check.names=FALSE)
 assign(paste(filenames[x]), data.in)
}
combined.df <-NULL
combined.df <- data.frame(row.names=row.names(ES09))
combined.df$ES09max <- apply(ES09[135:175], 1, max)
combined.df$ES10max <- apply(ES10[135:175], 1, max)
combined.df$ES11max <- apply(ES11[135:175], 1, max)
combined.df$ES12max <- apply(ES12[135:175], 1, max)
combined.df$ES13max <- apply(ES13[135:175], 1, max)
combined.df$ES14max <- apply(ES14[135:175], 1, max)
combined.df$ES15max <- apply(ES15[135:175], 1, max)
combined.df$ES16max <- apply(ES16[135:175], 1, max)

combined.df$ES09mean <- apply(ES09[135:175], 1, mean)
combined.df$ES10mean <- apply(ES10[135:175], 1, mean)
combined.df$ES11mean <- apply(ES11[135:175], 1, mean)
combined.df$ES12mean <- apply(ES12[135:175], 1, mean)
combined.df$ES13mean <- apply(ES13[135:175], 1, mean)
combined.df$ES14mean <- apply(ES14[135:175], 1, mean)
combined.df$ES15mean <- apply(ES15[135:175], 1, mean)
combined.df$ES16mean <- apply(ES16[135:175], 1, mean)
#temp <- ES09[135:175]
#################
## NORMALISING ##
#################

factors <- data.frame(means = colMeans(combined.df))
factors$variances <- colVars(as.matrix(combined.df))

ES09norm <- data.frame(apply(ES09, 2, function(x) (x - factors[ "ES09mean", "means" ]) / factors[ "ES09mean", "variances" ]))
ES10norm <- data.frame(apply(ES10, 2, function(x) (x - factors[ "ES10mean", "means" ]) / factors[ "ES10mean", "variances" ]))
ES11norm <- data.frame(apply(ES11, 2, function(x) (x - factors[ "ES11mean", "means" ]) / factors[ "ES11mean", "variances" ]))
ES12norm <- data.frame(apply(ES12, 2, function(x) (x - factors[ "ES12mean", "means" ]) / factors[ "ES12mean", "variances" ]))
ES13norm <- data.frame(apply(ES13, 2, function(x) (x - factors[ "ES13mean", "means" ]) / factors[ "ES13mean", "variances" ]))
ES14norm <- data.frame(apply(ES14, 2, function(x) (x - factors[ "ES14mean", "means" ]) / factors[ "ES14mean", "variances" ]))
ES15norm <- data.frame(apply(ES15, 2, function(x) (x - factors[ "ES15mean", "means" ]) / factors[ "ES15mean", "variances" ]))
ES16norm <- data.frame(apply(ES16, 2, function(x) (x - factors[ "ES16mean", "means" ]) / factors[ "ES16mean", "variances" ]))
colnames(ES09norm) <- colnames(ES09)
colnames(ES10norm) <- colnames(ES10)
colnames(ES11norm) <- colnames(ES11)
colnames(ES12norm) <- colnames(ES12)
colnames(ES13norm) <- colnames(ES13)
colnames(ES14norm) <- colnames(ES14)
colnames(ES15norm) <- colnames(ES15)
colnames(ES16norm) <- colnames(ES16)

combined.df$ES09normmax <- apply(ES09norm[135:175], 1, max)
combined.df$ES10normmax <- apply(ES10norm[135:175], 1, max)
combined.df$ES11normmax <- apply(ES11norm[135:175], 1, max)
combined.df$ES12normmax <- apply(ES12norm[135:175], 1, max)
combined.df$ES13normmax <- apply(ES13norm[135:175], 1, max)
combined.df$ES14normmax <- apply(ES14norm[135:175], 1, max)
combined.df$ES15normmax <- apply(ES15norm[135:175], 1, max)
combined.df$ES16normmax <- apply(ES16norm[135:175], 1, max)

combined.df$ES09normmean <- apply(ES09norm[135:175], 1, mean)
combined.df$ES10normmean <- apply(ES10norm[135:175], 1, mean)
combined.df$ES11normmean <- apply(ES11norm[135:175], 1, mean)
combined.df$ES12normmean <- apply(ES12norm[135:175], 1, mean)
combined.df$ES13normmean <- apply(ES13norm[135:175], 1, mean)
combined.df$ES14normmean <- apply(ES14norm[135:175], 1, mean)
combined.df$ES15normmean <- apply(ES15norm[135:175], 1, mean)
combined.df$ES16normmean <- apply(ES16norm[135:175], 1, mean)

samples = c("ES09","ES10","ES11","ES12","ES13","ES14","ES15","ES16","ES09","ES10","ES11","ES12","ES13","ES14","ES15","ES16","ES09","ES10","ES11","ES12","ES13","ES14","ES15","ES16","ES09","ES10","ES11","ES12","ES13","ES14","ES15","ES16")
type = c("max","max","max","max","max","max","max","max","mean","mean","mean","mean","mean","mean","mean","mean","max","max","max","max","max","max","max","max","mean","mean","mean","mean","mean","mean","mean","mean")
corr = c("raw","raw","raw","raw","raw","raw","raw","raw","raw","raw","raw","raw","raw","raw","raw","raw","norm","norm","norm","norm","norm","norm","norm","norm","norm","norm","norm","norm","norm","norm","norm","norm")

#combined.df2 <- data.frame(t(rbind("samples" = samples, "type" = type,"corr"=corr, combined.df)))


#########################
## Normalisation range ##
#########################
#refgenes.df <- unique(grep(paste(toMatch,collapse="|"), combined.df$gene, value=TRUE))

combined.std.df <- decostand(combined.df[1:32], 'standardize', MARGIN=2)
combined.dec.df <- decostand(combined.std.df, 'range', MARGIN=2)

################################
## RNAseq data (HTseq output) ##
################################

combined.df$gene <- row.names(combined.df)
combined.dec.df$gene <-  row.names(combined.dec.df)

# Significantly different RNA expression genese (FDR filtered)
#HTseq <- read.csv("/home/sbi6dap/Projects/ALD/RNAseq/ARA11/HTseq/alltags_edgeR_lt0.5FDR.csv")
# All expression comparisons
HTseq <- read.csv("/home/sbi6dap/Projects/ALD/RNAseq/ARA11/HTseq/alltags_edgeR.csv")
RNAraw <- read.table("/home/sbi6dap/Projects/ALD/RNAseq/ARA11/HTseq/all.count", header=TRUE)

DRcombined.dec.df <- merge(DRcombined.dec.df, HTseq, by = "gene", incomparables = NA)
DRcombined.dec.df <- merge(DRcombined.dec.df, RNAraw, by = "gene", incomparables = NA)

DRcombined.df <- merge(combined.df, HTseq, by = "gene", incomparables = NA)
DRcombined.df <- merge(combined.df, RNAraw, by = "gene", incomparables = NA)

#############
## Testing ##
#############
## FOLD CHANGE
# raw calcs
DRcombined.dec.df$raw.max.fold<- apply(DRcombined.dec.df[,2:45], 1, function (x) (mean(x[c(1,2,5,6)]) / mean(x[c(3,4,7,8)])))

# normalised calcs
DRcombined.dec.df$norm.max.fold <- apply(DRcombined.dec.df[,2:45], 1, function (x) (mean(x[c(17,18,21,22)]) / mean(x[c(19,20,23,24)])))
summary(DRcombined.dec.df$norm.max.fold)

# extract <0.5 & >2
DRcombined.dec.fold <- DRcombined.dec.df[ which(DRcombined.dec.df$norm.max.fold < 0.5 | DRcombined.dec.df$norm.max.fold > 2), ]

## T-TESTING
# raw calcs
raw.t.result <- apply(DRcombined.dec.df[,2:45], 1, function (x) t.test(x[c(1,2,5,6)],x[c(3,4,7,8)],paired=TRUE))
DRcombined.dec.df$rawmax.tpval <- unlist(lapply(raw.t.result, function(x) x$p.value))
DRcombined.dec.df$rawmax.fdr <- p.adjust(DRcombined.dec.df$rawmax.tpval, method = "fdr")

# normalised calcs
normmax.t.result <- apply(DRcombined.dec.df[,2:45], 1, function (x) t.test(x[c(17,18,21,22)],x[c(19,20,23,24)],paired=TRUE))
DRcombined.dec.df$normmax.tpval <- unlist(lapply(normmax.t.result, function(x) x$p.value))
DRcombined.dec.df$normmax.fdr <- p.adjust(DRcombined.dec.df$normmax.tpval, method = "fdr")
summary(DRcombined.dec.df$normmax.fdr)

# Annotate based on RNA expression
DRcombined.sigdiff <- DRcombined.dec.df[ which(DRcombined.dec.df$normmax.tpval < 0.05), ]
combined.df$logFC.quart <- c("Low","LQ","UQ","High")[ findInterval(combined.df$logFC, c(-Inf, 0.4725, 0.752, 1.3040, Inf))]

#############
### Charting
#############
library(ggplot2)
library(reshape2)
library(scales)

p <- ggplot(data=DRcombined.dec.df, aes(alpha=0.5), na.rm=TRUE)
p +
  geom_point(aes(x=log(ES09normmax + 1), y=log(ES3)), colour="red") + 
  geom_point(aes(x=log(ES10normmax + 1), y=log(ES3)), colour="red") + 
  geom_point(aes(x=log(ES11normmax + 1), y=log(ES4))) +
  geom_point(aes(x=log(ES12normmax + 1), y=log(ES4))) +
  geom_point(aes(x=log(ES13normmax + 1), y=log(ES5)), colour="red") +
  geom_point(aes(x=log(ES14normmax + 1), y=log(ES5)), colour="red") +
  geom_point(aes(x=log(ES15normmax + 1), y=log(ES6))) +
  geom_point(aes(x=log(ES16normmax + 1), y=log(ES6))) 
  coord_cartesian(xlim = c(0,0.01)) +
  labs(x="Chromatin occupancy summit (lt120, TSS +/- 250bp)",y="RNA HTseq count")
  #stat_summary(fun.y = mean, geom="line", aes(x=log(ES09max), y=log(ES3)))

p +
  geom_smooth(aes(x=log(ES09normmax), y=log(ES3 + 1)), se=FALSE, col="red", lty=2) +
  geom_smooth(aes(x=log(ES10normmax), y=log(ES3 + 1)), se=FALSE, col="blue", lty=2) + 
  geom_smooth(aes(x=log(ES11normmax), y=log(ES4 + 1)), se=FALSE, col="red", lty=1) +
  geom_smooth(aes(x=log(ES12normmax), y=log(ES4 + 1)), se=FALSE, col="blue", lty=1) +
  geom_smooth(aes(x=log(ES13normmax), y=log(ES5 + 1)), se=FALSE, col="red", lty=2) +
  geom_smooth(aes(x=log(ES14normmax), y=log(ES5 + 1)), se=FALSE, col="blue", lty=2) +
  geom_smooth(aes(x=log(ES15normmax), y=log(ES6 + 1)), se=FALSE, col="red", lty=1) +
  geom_smooth(aes(x=log(ES16normmax), y=log(ES6 + 1)), se=FALSE, col="blue", lty=1) +
  #coord_cartesian(xlim = c(0,0.125)) +
  labs(x="raw chromatin occupancy summit (lt120, CSS +/- 250bp)",y="RNA HTseq count")

############
# Boxplots #
############

combined.melt <-melt(combined.df, id=c("samples","type","corr","logFC.quart"))
combined.melt.raw <- subset(combined.melt, grepl("raw", corr))
combined.melt.norm <- subset(combined.melt, grepl("norm", corr))

p1 <- ggplot(data=subset.melt)
p2 <- ggplot(data=combined.melt.norm)

mp1 <- p1 +
  geom_boxplot(aes(x=samples, y=log(as.numeric(value)), fill=type), outlier.shape=NA) +
  coord_cartesian(ylim = c(-0.5,4)) + 
  labs(title = "raw")
  
mp2 <- p2 +
  geom_boxplot(aes(x=samples, y=log(as.numeric(value)), fill=type), outlier.shape=NA) +
  coord_cartesian(ylim = c(-11,-1)) +
  labs(title = "normalised")

multiplot(p1,p2)

subset.df <- data.frame(combined.df[18:25])
subset.df$logFC.quart <- combined.df$logFC.quart
subset.melt <-melt(subset.df, id=c("logFC.quart"))

p3 <- ggplot(data=subset.melt, aes(x=variable, y=log(as.numeric(value)), fill=logFC.quart))
p3 +
  geom_boxplot(outlier.shape=NA) 
  coord_cartesian(ylim = c(-9.5,-1)) +
  labs(title = "normalised")

###############
## misc plot ##
###############
tmp <- data.frame(t(rbind("raw" = ES11[120,], "norm" = ES11norm[120,])))
tmp2 <- tmp[135:175, ]
tmp.dec <- decostand(tmp2, 'range', MARGIN=2)
tmp.dec$pos <- row.names(tmp[135:175, ])


tmp.melt <- melt(tmp.dec)

p <- ggplot(data=tmp.melt, aes(x=as.numeric(pos), y=value), colour=variable)
p +
  stat_smooth(method="loess", span=0.1, se=FALSE) +
  scale_x_continuous(breaks = pretty_breaks(n=12))

library(ggrepel)

q <- ggplot(data=DRcombined.dec.df, aes(x=log(norm.max.fold), y=logFC, alpha=0.2))
q + 
  geom_point(data = subset(DRcombined.dec.df, FDR <0.05), colour="red") +
  geom_point(data = subset(DRcombined.dec.df, FDR >0.05), colour="black") +
  geom_smooth(data = subset(DRcombined.dec.df, FDR <0.05), se=TRUE, col="blue", lty=1)
  geom_smooth(se=TRUE, col="blue", lty=2) 
  geom_text_repel(data = subset(DRcombined.dec.fold, abs(logFC * log(norm.max.fold)) > 2), aes(label = gene)) 

