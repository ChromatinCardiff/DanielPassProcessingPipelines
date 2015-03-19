install.packages("ggplot2")
library(ggplot2)

x <- read.table("~/Projects/ACS/original_data/2012_2013Prelim_Ang_Marg/Exp1/split_chr/atfull_3peat.txt", header=FALSE)
colnames(x) <-c("chr", "nuclAbun", "pos", "strand", "rnaCount", "shape", "type", "gene")
#options(scipen=5)

#x.not0 <- subset(x, !grepl(0, val))

x.positive <- subset(x, grepl("\\+", strand))
x.utr <- subset(x.positive, grepl("5\\'utr", type))

x.sort <- x.utr[order(x.utr$pos),]

x.agg <-aggregate(x.sort, list(Region = x.sort$pos), mean)

x.narrow <- subset(x.sort, grepl("NarrowPeak", shape))
x.narrow_sub <- x.narrow[sample(nrow(x.narrow), 54477), ]
x.n.agg <-aggregate(x.narrow_sub, list(Region = x.narrow_sub$pos), mean)

x.broad <- subset(x.sort, grepl("BroadWithPeak", shape))
x.broad_sub <- x.broad[sample(nrow(x.broad), 54477), ]
x.b.agg <-aggregate(x.broad_sub, list(Region = x.broad_sub$pos), mean)

x.weak <- subset(x.sort, grepl("WeakPeak", shape))
x.weak_sub <- x.weak[sample(nrow(x.weak), 54477), ]
x.w.agg <-aggregate(x.weak_sub, list(Region = x.weak_sub$pos), mean)

#line plot overlays
#p <- ggplot()
#p + 
#  geom_ribbon(data=x.agg, aes(x=pos, ymin=0, ymax=nuclAbun)) +
#  geom_line(data=x.n.agg, aes(x=pos, y=nuclAbun), colour="blue") +
#  geom_line(data=x.b.agg, aes(x=pos, y=nuclAbun), colour="red") +
#  geom_line(data=x.w.agg, aes(x=pos, y=nuclAbun), colour="green")
#dev.off()

##facet plots
png("./At_full_normalised_5'UTR_region.png", width = 1920, height = 960, units = "px")
p <- ggplot() + theme(panel.grid.major = element_line(size = .5, color = "grey")) + coord_cartesian(ylim = c(0, 25), xlim = c(-500, 2000)) + geom_vline(xintercept = 0, colour="blue", size = 2, linetype="dotted")
p1 <- p + geom_ribbon(data=x.agg, aes(x=pos, ymin=0, ymax=nuclAbun)) + labs(title = "At genome - All Data", x="Distance from 5'UTR modal peak", y="Nucleosome abundance count")
p2 <- p + geom_ribbon(data=x.n.agg, aes(x=pos, ymin=0, ymax=nuclAbun)) + labs(title = "At genome - Narrow 5'UTR Sites", x="Distance from 5'UTR modal peak", y="Nucleosome abundance count")
p3 <- p + geom_ribbon(data=x.b.agg, aes(x=pos, ymin=0, ymax=nuclAbun)) + labs(title = "At genome - Broad 5'UTR Sites", x="Distance from 5'UTR modal peak", y="Nucleosome abundance count")
p4 <- p + geom_ribbon(data=x.w.agg, aes(x=pos, ymin=0, ymax=nuclAbun)) + labs(title = "At genome - Weak 5'UTR Sites", x="Distance from 5'UTR modal peak", y="Nucleosome abundance count")

multiplot(p1,p2,p3,p4, cols=2)
dev.off()

x.gene <- subset(x, grepl("AT5G36330", gene))
