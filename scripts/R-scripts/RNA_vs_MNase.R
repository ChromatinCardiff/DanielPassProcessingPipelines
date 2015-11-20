library(ggplot2)
library(reshape2)
library(scales)

setwd("~/Projects/ALD/RNAseq/ARA11/HTseq")

x <- read.table("samples.txt", header=TRUE)

x$LL <- (x$ES09 + x$ES13) / 2
x$DL <- (x$ES11 + x$ES15) / 2
x$LH <- (x$ES10 + x$ES14) / 2
x$DH <- (x$ES12 + x$ES16) / 2
head(x,10)

y <- x[c("Gene","LL","DL","LH","DH")]


y$Lowdiff <- x$LL - x$DL
y$Highdiff <- x$LH - x$DH

y.high <- y.high[with(y.high, order(diff)), ]

write.table(y, "d.all_summary.csv", sep=",")

x <- read.csv("/home/sbi6dap/Projects/ALD/RNAseq/TAIR10/edgeR_analysis/MNase-RNAseq.csv", header=TRUE)

x$logHighdiff <- log10(abs(x$Highdiff) +1)
x$logLowdiff <- log10(abs(x$Lowdiff) +1)


x.melt <-melt(x, id=c("Gene"))

p <- ggplot(data=x, aes(x=abs(logFC), y=logLowdiff))
p + geom_point() +
  geom_text(aes(label=ifelse(logLowdiff>1 & abs(logFC)>1,as.character(Gene),'')),hjust=0,just=0)

  stat_smooth(se=FALSE) 

  coord_cartesian(ylim=c(0,50))
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_vline(x=0, colour="black", lty=2) +
  geom_vline(x=168, colour="black", lty=2) +
  scale_colour_brewer(palette="Set1") +
  facet_wrap(~ Treatment, ncol=1)