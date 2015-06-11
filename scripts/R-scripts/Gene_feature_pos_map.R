library(ggplot2)
library(data.table)
library(scales)
library(reshape2)
library(plyr)
install.packages("vegan")
library(vegan)

# Read in
x <- read.table("~/Projects/ALD/MNase-seq/old_vs_new/combined_TSS.txt", header=TRUE, sep=",")
summary(x)

x1 <- read.table("~/Projects/ACS/analysis/dpos/profile_TSS_heatmap/TSS_full_uniq.txt", header=TRUE)
rownames(x1) <-x1$name
x1[,1] <- NULL
x2<-colMeans(x1, na.rm=TRUE)
x2$pos <-as.numeric(colnames(y1))
x.melt <-melt(x2)
x.melt2 <-na.omit(x.melt)

# Read extra datasets
y1 <- read.table("~/Projects/ACS/analysis/dpos/profile_TSS_heatmap/TSS_full_uniq_masked.txt", header=TRUE)
rownames(y1) <-y1$name
y1[,1] <- NULL
y2<-colMeans(y1, na.rm=TRUE)
y2$pos <-as.numeric(colnames(y1))
y.melt <-melt(y2)
y.melt2 <-na.omit(y.melt)


# Processing
rownames(x) <-x$pos
x.dec <- decostand(x, 'standardize', MARGIN=2)
x.dec$pos <- as.numeric(rownames(x))
x.melt <- melt(x.dec, id=c("pos"))

# Individual biorep plots
p <- ggplot(data=x.dec)
p1 <- p + 
  geom_line(aes(x=pos, y=ES09_150)) +
  geom_line(aes(x=pos, y=ES13_150)) +
  geom_line(aes(x=pos, y=Marg)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_vline(x=0, colour="blue") +
  labs(title = "Light Grown, Low MNase", x = "Position from TSS", y = "Normalised abundance (Hellinger)")

p2 <- p + 
  geom_line(aes(x=pos, y=ES10_150)) +
  geom_line(aes(x=pos, y=ES14_150)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_vline(x=0, colour="blue") +
  labs(title = "Light Grown, High MNase", x = "Position from TSS", y = "Normalised abundance (Hellinger)")

p3 <- p + 
  geom_line(aes(x=pos, y=ES11_150)) +
  geom_line(aes(x=pos, y=ES15_150)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_vline(x=0, colour="blue") +
  labs(title = "Dark Grown, Low MNase", x = "Position from TSS", y = "Normalised abundance (Hellinger)")

p4 <- p + 
  geom_line(aes(x=pos, y=ES12_150)) +
  geom_line(aes(x=pos, y=ES16_150)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_vline(x=0, colour="blue") +
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
  geom_vline(x=0, colour="blue") +
  labs(title = "Low MNase, Light vs Dark", x = "Position from TSS", y = "Normalised abundance (f=standardize)")

# High MNase Light vs Dark
p + 
  geom_ribbon(aes(x=pos, ymin=ES10_150, ymax=ES14_150, alpha=0.9, fill="Light|High")) +
  geom_ribbon(aes(x=pos, ymin=ES12_150, ymax=ES16_150, alpha=0.9, fill="Dark|High")) +
  geom_line(aes(x=pos, y=ES10_150)) +
  geom_line(aes(x=pos, y=ES14_150)) +
  geom_line(aes(x=pos, y=ES12_150)) +
  geom_line(aes(x=pos, y=ES16_150)) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_vline(x=0, colour="blue") +
  scale_colour_brewer(palette="Dark2") +
  labs(title = "HIGH MNase, Light vs Dark", x = "Position from TSS", y = "Normalised abundance (f=standardize)")

# All datasets
p2 <-ggplot(data=x.melt)
p2 +
  geom_line(aes(x=pos, y=value, colour=variable, )) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  scale_colour_brewer(palette="Dark2") +
  geom_vline(x=0, colour="blue")

p + 
  geom_line(aes(x=pos, y=ES09_150), colour="Red") +      # Low
  geom_line(aes(x=pos, y=ES13_150), colour="Red") +     # Low
  geom_line(aes(x=pos, y=ES10_150), colour="Green") +     # High
  geom_line(aes(x=pos, y=ES14_150), colour="Green") +     # High
  geom_line(aes(x=pos, y=Marg_150), colour="Black", size=2) +
  geom_line(aes(x=pos, y=AM1_150), colour="Darkblue", size=2) +
  geom_line(aes(x=pos, y=AM2_150), colour="Darkblue", size=2) +
  geom_line(aes(x=pos, y=AM3_150), colour="Darkblue", size=2) +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_vline(x=0, colour="blue") +
  labs(title = "Light Grown", x = "Position from TSS", y = "Normalised abundance (Hellinger)")


#### MISC
y.melt <- melt(y2)
y.melt
p <- ggplot(data=y.melt, aes(x=Var1, y=value))
p + 
  geom_line() +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  geom_vline(x=0, colour="blue") +
  labs(title = "'Marg'", x = "Position from TSS")


## Single gene
