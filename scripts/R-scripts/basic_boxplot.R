x <- read.table("~/temp/minus1.tsv", header=TRUE, sep="\t", row.names=1)

x.melt <- melt(x)
x.melt

minus1p <- t.test(x$minus1~x$digest)$p.value
plus1p <- t.test(x$plus1~x$digest)$p.value

ggplot(x.melt, aes(x=variable, y=value, fill=digest)) +
  geom_boxplot() +
  scale_fill_brewer(palette="Dark2") +
  # Pvalues
  # annotate("text", x=1, y=1.225, label= paste("pvalue =", round(minus1p, digits=3))) +
  #annotate("text", x=2, y=1.445, label= paste("pvalue =",  round(plus1p, digits=3))) +
  # Alternative
  annotate("text", x=1, y=1.225, label= paste("*"), size=10) +
  annotate("text", x=2, y=1.445, label= paste("*"), size=10) +
  # bars
  #1st
  annotate("segment", x=0.81, y=1.2, xend=1.19, yend=1.2) +
  annotate("segment", x=0.81, y=1.2, xend=0.81, yend=1.18) +
  annotate("segment", x=1.19, y=1.2, xend=1.19, yend=1.18) +
  #2nd
  annotate("segment", x=1.81, y=1.44, xend=2.19, yend=1.44) +
  annotate("segment", x=1.81, y=1.44, xend=1.81, yend=1.42) +
  annotate("segment", x=2.19, y=1.44, xend=2.19, yend=1.42) 
