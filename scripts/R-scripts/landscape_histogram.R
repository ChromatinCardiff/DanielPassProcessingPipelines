library(ggplot2)

x <- read.table("~/Projects/ACS/original_data/2012_2013Prelim_Ang_Marg/Exp1/Athal_4_Part150_AB_10.sgr", header=FALSE)

summary(x)

x.ch4 <- subset(x, grepl("chr4", V1))

ggplot(x.ch4, aes(y=log10(V3), x=V2)) +
  geom_histogram(=0.1, stat="identity") 
  facet_wrap(~V1, ncol = 1)
