x.dec$pos <- rownames(x.dec)
x.deccut <- x.dec[100:150,]
x.melt <-melt(x.dec, "pos")

p <-ggplot(data=x.melt, aes(x=as.numeric(as.character(pos)), y=value, colour=variable))
p +
  stat_smooth(method="loess", span=0.1, se=FALSE) 
  #geom_point() 
  scale_x_continuous(breaks = pretty_breaks(n=12)) +
  #scale_colour_brewer(palette="Paired") +
  geom_vline(x=0, colour="blue", lty=2) +
  #scale_y_continuous(limits=c(-3.5, 3.5)) +
  theme(legend.position = "bottom") +
  labs(title = "TSS") 
  facet_wrap(~ Treatment, ncol=1)
TSS.plot


####


tmp <- read.table("/home/sbi6dap/Projects/ALD/totalcoverage/RAWDATA/full_totalcov_num.txt", header=FALSE, sep="\t")

library(assertthat)

winsorize <-
  function(x, q=0.05)
  {
    assert_that(is.numeric(x))
    assert_that(is.number(q), q>=0, q<=1)
    
    lohi <- quantile(x, c(q, 1-q), na.rm=TRUE)
    if(diff(lohi) < 0) lohi <- rev(lohi)
    
    #x[!is.na(x) & x < lohi[1]] <- lohi[1]
    x[!is.na(x) & x > lohi[2]] <- lohi[2]
    x
  }

tmp.win <-tmp
tmp.win[, -1] <- sapply(tmp.win[,-1], winsorize)
tmp.dec <- data.frame(decostand(tmp.win, 'normalize', MARGIN=2))

tmp.dec$abun <- ((tmp.dec$V1 + tmp.dec$V3 + tmp.dec$V5 + tmp.dec$V7) / 4) - ((tmp.dec$V2 + tmp.dec$V4 + tmp.dec$V6 + tmp.dec$V8) / 4) # Low - High
write.table(tmp.dec$abun, "/home/sbi6dap/Projects/ALD/totalcoverage/low-high_digest.wig", sep="\t")

