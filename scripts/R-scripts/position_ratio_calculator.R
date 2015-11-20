library(ggplot2)
library(reshape2)
library(scales)

setwd("~/Projects/ALD/totalcoverage/mean_CSS_heatmap")

ES09 <- read.table("ES09_totalcov.xls.cut", header=TRUE, sep="\t", row.names=1)
ES10 <- read.table("ES10_totalcov.xls.cut", header=TRUE, sep="\t", row.names=1)
ES11 <- read.table("ES11_totalcov.xls.cut", header=TRUE, sep="\t", row.names=1)
ES12 <- read.table("ES12_totalcov.xls.cut", header=TRUE, sep="\t", row.names=1)
ES13 <- read.table("ES13_totalcov.xls.cut", header=TRUE, sep="\t", row.names=1)
ES14 <- read.table("ES14_totalcov.xls.cut", header=TRUE, sep="\t", row.names=1)
ES15 <- read.table("ES15_totalcov.xls.cut", header=TRUE, sep="\t", row.names=1)
ES16 <- read.table("ES16_totalcov.xls.cut", header=TRUE, sep="\t", row.names=1)

# Annotations
samples <- list(ES09,ES10,ES11,ES12,ES13,ES14,ES15,ES16)
samplenames <- list("ES09","ES10","ES11","ES12","ES13","ES14","ES15","ES16")
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

ratio.df <- NULL
for(i in samples){
  ratio.df <- cbind(ratio.df, CalculateRatio(i,140,160,150))
}
colnames(ratio.df) <- c(samplenames)

res <- Reduce(function(a,b){
  ans <- merge(a,b,by="row.names",incomparables = NA)
  row.names(ans) <- ans[,"Row.names"]
  ans[,!names(ans) %in% "Row.names"]
}, list(all.count, HTseq, ratio.df))

write.table(res, "/home/sbi6dap/Projects/ALD/totalcoverage/LDgenes_totalcov_vs_exp.txt", sep="\t")

resFilter <- res[abs(res$FDR) < 0.05,]
resFilterInv <- res[abs(res$FDR) > 0.05,]

tmp.melt <- melt(res)
#######################################

p <-ggplot(data=tmp.melt)
p +
  stat_smooth(aes(x=as.numeric(as.character(variable)), y=value, method="loess", span=0.01, se=FALSE) + 
  #geom_line(aes(x=as.numeric(as.character(variable)), y=value)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  #scale_colour_brewer(palette="Paired") +
  geom_vline(x=0, colour="blue", lty=2) +
  #scale_y_continuous(limits=c(-3.5, 3.5)) +
  theme(legend.position = "bottom") +
  labs(title = "gene")

plot(data=resFilter, log(ES09) ~ log(ES3))
plot(data=resFilterInv, log(ES09) ~ log(ES3))
