library(ggplot2)
library(reshape2)
library(scales)
library(vegan)


setwd("~/Projects/ALD/z/mean_CSS_heatmap")

ES09 <- read.table("ES09_totalcov.xls.cut", header=TRUE, sep="\t", row.names=1, check.names = FALSE)
ES10 <- read.table("ES10_totalcov.xls.cut", header=TRUE, sep="\t", row.names=1)
ES11 <- read.table("ES11_totalcov.xls.cut", header=TRUE, sep="\t", row.names=1)
ES12 <- read.table("ES12_totalcov.xls.cut", header=TRUE, sep="\t", row.names=1)
ES13 <- read.table("ES13_totalcov.xls.cut", header=TRUE, sep="\t", row.names=1)
ES14 <- read.table("ES14_totalcov.xls.cut", header=TRUE, sep="\t", row.names=1)
ES15 <- read.table("ES15_totalcov.xls.cut", header=TRUE, sep="\t", row.names=1)
ES16 <- read.table("ES16_totalcov.xls.cut", header=TRUE, sep="\t", row.names=1)

# Annotations
samplesHigh <- list(ES10,ES12,ES14,ES16)
samplesLow <- list(ES09,ES11,ES13,ES15)
samplesNamesHigh <- list("ES10","ES12","ES14","ES16")
samplesNamesLow <- list("ES09","ES11","ES13","ES15")
sampleNames <- list("ES09","ES10","ES11","ES12","ES13","ES14","ES15","ES16")
Exposure <- c("Light","Light","Dark","Dark","Light","Light","Dark","Dark")
Treatment <- c("Low","High","Low","High","Low","High","Low","High")
datatype <- c("datatype","Chromatin","Chromatin","Chromatin","Chromatin","Chromatin","Chromatin","Chromatin","Chromatin","RNAexp","RNAexp","RNAexp","RNAexp","RNAexp","RNAexp","RNAexp","RNAexp","RNAexp","RNAexp","RNAexp","RNAexp")

# RNA count data
all.count <- read.table("~/Projects/ALD/RNAseq/ARA11/HTseq/all.count", header=TRUE, sep="\t", row.names=1)
HTseq <- read.table("~/Projects/ALD/RNAseq-annotation-results/alltags/alltags_edgeR.csv", header=TRUE, sep=",", row.names=1)


CalculateRatio <- function(Sample,Left,Right,Split)
{
  MeanLeft <- rowMeans(Sample[Left:Split])
  MeanRight <- rowMeans(Sample[Split:Right])
  
  ratio <- MeanLeft/MeanRight
  ratio
}

Lowratio.df <- NULL
Lowmean.df <- NULL

for(i in samplesLow){
  Lowratio.df <- data.frame(cbind(Lowratio.df, CalculateRatio(i,110,190,150)))
  Lowmean.df <- data.frame(cbind(Lowmean.df, rowMeans(i[110:150])))
}

colnames(Lowratio.df) <- c(samplesNamesLow)
Lowratio.df$Light <- ((Lowratio.df$ES09 + Lowratio.df$ES13)/2)
Lowratio.df$Dark <- ((Lowratio.df$ES11 + Lowratio.df$ES15)/2)
Lowratio.df$diff <- Lowratio.df$Light - Lowratio.df$Dark

resratios <- Reduce(function(a,b){
  ans <- merge(a,b,by="row.names",incomparables = NA)
  row.names(ans) <- ans[,"Row.names"]
  ans[,!names(ans) %in% "Row.names"]
}, list(HTseq, Lowratio.df))                             # all.count,


#write.table(downreg.list, "/home/sbi6dap/Projects/ALD/totalcoverage/downreg_full.txt", sep="\t")

resratios$gene <-row.names(resratios)

resratiosFilter <- resratios[abs(resratios$FDR) < 0.05,]
resratiosFilter2 <- resratiosFilter[abs(resratiosFilter$logFC) > 2,]
resFilterInv <- res[abs(res$FDR) > 0.05,]

tmp.melt <- melt(resratiosFilter, id=c("logFC","logCPM","PValue","FDR","gene"))
tmp.melt.sub <- subset(tmp.melt, grepl("diff", variable))
tmp.melt.filter <- tmp.melt[abs(tmp.melt$value) > 1,]
#######################################

p <-ggplot(data=tmp.melt.sub, aes(x=value, y=abs(logFC), alpha=0.5))  # colour=FDR,
p +
  labs(title = "Difference in utr occupancy upstream of CDS vs gene foldchange. Only genes >0.05 FDR represented ", x="Ratio between 400bp upstream of CDS to 400 downstream, light - dark", y="abs RNA logFC") +
  #coord_cartesian(xlim =c(-10, 100), ylim =c(0, 1)) +
  xlim(-1,1) +
  #scale_colour_brewer(palette="Paired") +
  geom_point() +
  geom_text(aes(label=ifelse(value>0.5 & abs(logFC)>1,as.character(gene),'')),hjust=0,just=0) +
  geom_text(aes(label=ifelse(value>0.25 & abs(logFC)>3,as.character(gene),'')),hjust=0,just=0) +
  geom_text(aes(label=ifelse(value< (-0.5) & abs(logFC)>1,as.character(gene),'')),hjust=1,just=0) +
  geom_text(aes(label=ifelse(value< (-0.25) & abs(logFC)>3,as.character(gene),'')),hjust=1,just=0)

  stat_smooth(aes(group=variable), method = loess, se=TRUE)

  geom_point(aes(x=Light, y=logFC), colour="red") +
  geom_point(aes(x=Dark, y=logFC), colour="blue") 
  geom_point(aes(x=ES11, y=logFC), colour="red") +
  geom_point(aes(x=ES12, y=logFC), colour="blue") +
  geom_point(aes(x=ES13, y=logFC), colour="red") +
  geom_point(aes(x=ES14, y=logFC), colour="blue") +
  geom_point(aes(x=ES15, y=logFC), colour="red") +
  geom_point(aes(x=ES16, y=logFC), colour="blue") +

  
## Normal chart
upreg.list <- resratios[resratios$diff > 0.5,]
uplist <-upreg.list$gene
downreg.list <- resratios[resratios$diff < -0.5,]
downlist <-downreg.list$gene

ES09$gene <- row.names(ES09)
posuptest <- ES09[match(upreg.list$gene, ES09$gene), ]
  #posup <- data.frame(merge(ES09, uplist, by = "row.names", incomparables = NA, check.names = FALSE))
  #posdown <- merge(ES09, downreg.list, by = "row.names", incomparables = NA)
posuptest$gene <- NULL
upchromatin.dec <- data.frame(decostand(posuptest, 'standardize', MARGIN=2))
upchromatin.means <- data.frame(mean=colMeans(upchromatin.dec))
upchromatin.means$pos <- row.names(upreg.list$gene)
up.melt <- melt(upchromatin.means, id=c("pos"))


############## BOXPLOTS ####################
allreg.melt <- melt(allreg.list, id=c("logFC","logCPM","PValue","FDR","gene","direction","diff"))
summary(allreg.list)
p <-ggplot(data=allreg.melt)
ylim1 = c(boxplot.stats(log(allreg.melt$logFC))$stats[3] - abs((boxplot.stats(log(allreg.melt$logFC))$stats[1]*1.1)),boxplot.stats(log(allreg.melt$value))$stats[3] + (boxplot.stats(log(allreg.melt$logFC))$stats[5]*1.1))
p +
  labs(title = "DE for genes with >0.05 5'UTR to Exon1 ratio difference", x="5'UTR to Exon occupancy change", y="Ratio ") +
  geom_boxplot(aes(x=as.factor(direction), y=logCPM, fill=variable), outlier.shape = NA) 
  #geom_boxplot(aes(x=as.factor(direction), y=log(Dark))) +
  coord_cartesian(ylim = ylim1)
  coord_cartesian(ylim = c(-0.5,0.5))
  #scale_colour_brewer(palette="Paired") +
  geom_vline(x=0, colour="blue", lty=2) 
  
  theme(legend.position = "bottom") 
  labs(title = "TSS") +
  facet_wrap(~ Treatment, ncol=1)

plot(data=resFilter, log(ES09) ~ log(ES3))
plot(data=resFilterInv, log(ES09) ~ log(ES3))

ES09up <- read.table("ES09_tc_up.xls", header=TRUE, sep="\t", row.names=1, check.names = FALSE)
ES09down <- read.table("ES09_tc_down.xls", header=TRUE, sep="\t", row.names=1, check.names = FALSE)

library(assertthat)

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

ES09up.win <- ES09up
ES09up.win[, -1] <- sapply(ES09up.win[,-1], winsorize)
ES09down.win <- ES09down
ES09down.win[, -1] <- sapply(ES09down.win[,-1], winsorize)

ES09up.dec <- data.frame(decostand(ES09up.win, 'standardize', MARGIN=1))
ES09down.dec <- data.frame(decostand(ES09down.win, 'standardize', MARGIN=1))

ES09up.means <- data.frame(meanup=colMeans(ES09up.dec))#, pos=as.numeric(as.character(colnames(ES09up.win))))
ES09down.means <- data.frame(meandown=colMeans(ES09down.dec))#, pos=as.numeric(as.character(colnames(ES09down.win))))

ES09.full <- cbind(ES09up.means,ES09down.means, pos=as.numeric(as.character(colnames(ES09down.win))))

ES09.melt <- melt(ES09.full, id=c("pos"))

p <-ggplot(data=ES09.melt, aes(x=pos, y=value), colour=variable)
p +
  #stat_smooth(method="loess", span=0.01, se=FALSE) + 
  geom_line(aes(colour=variable)) #+#, colour=Exposure)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  #scale_colour_brewer(palette="Paired") +
  geom_vline(x=0, colour="blue", lty=2) 
  #scale_y_continuous(limits=c(-3.5, 3.5)) +
  theme(legend.position = "bottom") 
  labs(title = "TSS") +
  facet_wrap(~ Treatment, ncol=1)


averageES09$xcol <- row.names(averageES09)
p+geom_point()
