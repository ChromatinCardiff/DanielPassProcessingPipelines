library(vegan)

setwd("/home/sbi6dap/Projects/ALD/totalcoverage/size_comparisons/ALL_TSS")
# Chromatin per-gene profile files
myFiles <- list.files(pattern="*DGE.xls-nodup")
filenames = c("ES09","ES11","ES13","ES15")

for(x in seq_along(myFiles)){
 data.in <- read.table(myFiles[x], header=TRUE, sep="\t", row.names=1)
 assign(paste(filenames[x]), data.in)
}

combined.df <- data.frame(row.names=row.names(ES09))
combined.df$ES09max <- apply(ES09[175:225], 1, max)
combined.df$ES11max <- apply(ES11[175:225], 1, max)
combined.df$ES13max <- apply(ES13[175:225], 1, max)
combined.df$ES15max <- apply(ES15[175:225], 1, max)

###########################
## Normalisation by gene ##
###########################
toMatch <- c("AT2G37620","AT3G53750","AT2G42090","AT5G09810","AT3G18780","AT5G59370","AT1G49240")
refgenes.df <- unique(grep(paste(toMatch,collapse="|"), combined.df$gene, value=TRUE))

refFactors <-c(1.7789398433,1.9361485337,4.4916923566,14.2508739881,9.2786804084,5.3402972284,11.5266550982,5.3699067185)

#combined.df$ES09norm <- combined.df$ES09max * 1
#combined.df$ES11norm <- combined.df$ES11max * 14.4916923566
#combined.df$ES13norm <- combined.df$ES13max * 14.2786804084
#combined.df$ES15norm <- combined.df$ES15max * 17.5266550982

combined.std.df <- decostand(combined.df, 'standardize', MARGIN=2)
combined.dec.df <- decostand(combined.std.df, 'range', MARGIN=2)

combined.dec.df$gene <- row.names(combined.dec.df)
combined.df$gene <- row.names(combined.df)


# RNAseq data (HTseq output)
HTseq <- read.csv("/home/sbi6dap/Projects/ALD/RNAseq/ARA11/HTseq/alltags_edgeR_lt0.5FDR.csv")
RNAraw <- read.table("/home/sbi6dap/Projects/ALD/RNAseq/ARA11/HTseq/all.count", header=TRUE)

combined.dec.df <- merge(combined.dec.df, HTseq, by = "gene", incomparables = NA)
combined.dec.df <- merge(combined.dec.df, RNAraw, by = "gene", incomparables = NA)

combined.df <- merge(combined.df, HTseq, by = "gene", incomparables = NA)
combined.df <- merge(combined.df, RNAraw, by = "gene", incomparables = NA)


#############
### Charting
#############
library(ggplot2)
library(ggrepel)
library(reshape2)
library(scales)

p <- ggplot(data=combined.dec.df, na.rm=TRUE)
p +
  #geom_point(aes(x=log(ES09max), y=log(ES3))) 
  geom_smooth(aes(x=ES09max, y=logFC), colour="red") +
  geom_smooth(aes(x=ES11max, y=logFC)) +
  geom_smooth(aes(x=ES13max, y=logFC), colour="red") +
  geom_smooth(aes(x=ES15max, y=logFC)) +
  labs(x="Chromatin occupancy summit (lt120, TSS +/- 250bp)",y="logFC")
  #stat_summary(fun.y = mean, geom="line", aes(x=log(ES09max), y=log(ES3)))

p +
  geom_smooth(aes(x=log(ES09norm + 1), y=log(ES09max +1)), colour="red") +
  geom_smooth(aes(x=log(ES11norm + 1), y=log(ES11max +1)), colour="blue") +
  geom_smooth(aes(x=log(ES13norm + 1), y=log(ES13max +1)), colour="green") +
  geom_smooth(aes(x=log(ES15norm + 1), y=log(ES15max +1)), colour="black") 

