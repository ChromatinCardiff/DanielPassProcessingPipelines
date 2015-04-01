library(ggplot2)
require(scales)
library(reshape2)
library(plyr)
library(vegan)

# Read in
x <- read.table("~/Projects/ALD/dpos_peaks/profile_TSS_heatmap/combined_TSS.txt", header=TRUE)
summary(x)

# Processing
rownames(x) <-x$pos
x.dec <- decostand(x, 'standardize', MARGIN=2)
x.dec$pos <- as.numeric(rownames(x))
x.melt <- melt(x, id=c("pos"))

# Individual biorep plots
p <- ggplot(data=x.dec)
p1 <- p + 
  geom_line(aes(x=pos, y=ES09_150)) +
  geom_line(aes(x=pos, y=ES13_150)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_abline(x=150, colour="blue") +
  labs(title = "Light Grown, Low MNase", x = "Position from TSS", y = "Normalised abundance (Hellinger)")

p2 <- p + 
  geom_line(aes(x=pos, y=ES10_150)) +
  geom_line(aes(x=pos, y=ES14_150)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_abline(x=150, colour="blue") +
  labs(title = "Light Grown, High MNase", x = "Position from TSS", y = "Normalised abundance (Hellinger)")

p3 <- p + 
  geom_line(aes(x=pos, y=ES11_150)) +
  geom_line(aes(x=pos, y=ES15_150)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_abline(x=150, colour="blue") +
  labs(title = "Dark Grown, Low MNase", x = "Position from TSS", y = "Normalised abundance (Hellinger)")

p4 <- p + 
  geom_line(aes(x=pos, y=ES12_150)) +
  geom_line(aes(x=pos, y=ES16_150)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_abline(x=150, colour="blue") +
  labs(title = "Dark Grown, High MNase", x = "Position from TSS", y = "Normalised abundance (Hellinger)")

multiplot(p1,p2,p3,p4, cols=2)


# Low MNase, Light vs Dark
p + 
  geom_ribbon(aes(x=pos, ymin=ES09_150, ymax=ES13_150, alpha=0.9, fill="Light|Low")) +
  geom_ribbon(aes(x=pos, ymin=ES11_150, ymax=ES15_150, alpha=0.9, fill="Dark|Low")) +
  geom_line(aes(x=pos, y=ES09_150)) +
  geom_line(aes(x=pos, y=ES13_150)) +
  geom_line(aes(x=pos, y=ES11_150)) +
  geom_line(aes(x=pos, y=ES15_150)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_abline(x=150, colour="blue") +
  labs(title = "Low MNase, Light vs Dark", x = "Position from TSS", y = "Normalised abundance (f=standardize)")

# High MNase Light vs Dark
p + 
  geom_ribbon(aes(x=pos, ymin=ES10_150, ymax=ES14_150, alpha=0.9, fill="Light|Low")) +
  geom_ribbon(aes(x=pos, ymin=ES12_150, ymax=ES16_150, alpha=0.9, fill="Dark|Low")) +
  geom_line(aes(x=pos, y=ES10_150)) +
  geom_line(aes(x=pos, y=ES14_150)) +
  geom_line(aes(x=pos, y=ES12_150)) +
  geom_line(aes(x=pos, y=ES16_150)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_abline(x=150, colour="blue") +
  scale_colour_brewer(palette="Dark2") +
  labs(title = "HIGH MNase, Light vs Dark", x = "Position from TSS", y = "Normalised abundance (f=standardize)")
