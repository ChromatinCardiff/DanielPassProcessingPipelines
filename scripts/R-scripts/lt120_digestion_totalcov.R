library(reshape2)
library(ggplot2)
library(vegan)


setwd("/home/sbi6dap/Projects/ALD/totalcoverage/size_comparisons/")
       
df <- read.table("prop_lt120_TSS.xls", header=TRUE, sep="\t")
ref <- data.frame(t(read.table("refgenes/actin3.csv", header=FALSE, sep="\t", row.names=1)))
ref$pos <- df$pos

df.std <- decostand(df, 'standardize', MARGIN=2)
df.std$pos <- df$pos

ref.std <- decostand(ref, 'standardize', MARGIN=2)
ref.std$pos <- ref$pos

df.melt <-melt(df, id="pos")
df.std.melt <-melt(df.std, id="pos")

ref.melt <-melt(ref, id="pos")
ref.std.melt <-melt(ref.std, id="pos")


  
p1 <- ggplot(data=df.melt) +
  stat_smooth(aes(x=pos, y=value, colour=variable), method="loess", span=0.1, se=FALSE) 
  facet_grid(. ~ variable) +
  theme(legend.position = "none") + labs(title = "Total lt120 raw")

p2 <- ggplot(data=df.std.melt) +
  stat_smooth(aes(x=pos, y=value, colour=variable), method="loess", span=0.1, se=FALSE) +
  facet_grid(. ~ variable) +
  theme(legend.position = "none") + labs(title = "Total lt120 standardized")

p3 <- ggplot(data=ref.melt) +
  stat_smooth(aes(x=pos, y=value, colour=variable), method="loess", span=0.1, se=FALSE) +
  facet_grid(. ~ variable) +
  theme(legend.position = "none") + labs(title = "actin3 lt120")

p4 <- ggplot(data=ref.std.melt) +
  stat_smooth(aes(x=pos, y=value, colour=variable), method="loess", span=0.1, se=FALSE) 
  facet_grid(. ~ variable) +
  theme(legend.position = "none") + labs(title = "actin3 lt120 standardised")


multiplot(p1,p2,p3,p4)

#Individual chart
p2 <- ggplot(data=df.std.melt) +
  stat_smooth(aes(x=pos, y=value, colour=variable), method="loess", span=0.1, se=FALSE) +
  #facet_grid(. ~ variable) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_vline(xintercept=0, colour="blue", linetype="longdash") +
  theme(legend.position = "none") + labs(title = "Total lt120 standardized")
p2
