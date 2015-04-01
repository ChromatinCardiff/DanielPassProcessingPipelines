library(ggplot2)
library(scales)
library(sitools)

x <- read.table("/home/sbi6dap/Projects/ALD/sgr_150/ES09_150.sgr", header=FALSE)
options(scipen=2)
summary(x)

x.sub10 <- subset(x, grepl("dpos10", V4), na.rm=TRUE)
x.sub30 <- subset(x, grepl("dpos30", V4), na.rm=TRUE)

ggplot(x, aes(x=as.factor(V3), y=V2)) +
  geom_bin2d(binwidth = c(1,301), aes(fill=as.factor(..count..))) +
  scale_fill_manual(values = c("#000000","#9E8400","#AB7402","#FAD000","#FBAB04","#FC8608","#FD610C","#FE3C10","#FF1714")) +
  scale_y_continuous(breaks = pretty_breaks(12)) +
  facet_wrap(~V1, nrow = 1)

ggplot(x, aes(x=V1, y=V2)) +
  stat_density(aes(ymax = ..density.., ymin= -..density..), geom = "ribbon", position = "identity") 
  coord_flip() +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  facet_wrap(~V1, nrow = 1)
  
ggplot(rbind(data.frame(x.sub10, group="gr10"), data.frame(x.sub30, group="gr30")), aes(x=V2)) +
  stat_density(aes(colour= group, alpha=0.5, ymax = ..count.., ymin= -..count..), trim=TRUE, geom = "ribbon", position = "identity") +
  coord_flip() +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  facet_wrap(~V1, nrow = 1)

ggplot(x, aes(x=V2)) +
  stat_density(aes(alpha=0.5, ymax = ..count.., ymin= -..count..), trim=TRUE, geom = "ribbon", position = "identity") +
  coord_flip() +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  facet_wrap(~V1, nrow = 1)
