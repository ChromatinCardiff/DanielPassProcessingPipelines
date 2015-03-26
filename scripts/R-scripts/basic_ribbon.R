library(ggplot2)
library(zoo)
require(scales)
library(plyr)

x <- read.table("~/Projects/ALD/ALD_histogram.txt", header=FALSE)
summary(x)

fn <- function(x) x/sum(x)

### Normalise column V3
x.nor <- ddply(x, "V1", transform, V3norm=fn(V3))

### SMOOTHING
x$av <- ave(x.nor$V3norm, x.nor$V1, FUN= function(x) rollmean(x, k=30, fill=NA))
x$lav <-log10(x$av)

summary(x)

p <- ggplot(data=x, aes(V2, lav, colour=V1))
p + 
  geom_line() +
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  scale_colour_brewer(palette="Dark2")
  
