install.packages("ggplot2")
library(ggplot2)

x <- read.table("~/Projects/ACS/original_data/2012_2013Prelim_Ang_Marg/Exp1/split_chr/chr1_23.sgr", header=FALSE)
colnames(x) <-c("chr", "nuclAbun", "pos", "strand", "rnaCount", "shape", "type", "gene")
#options(scipen=5)

#x.not0 <- subset(x, !grepl(0, val))

x.positive <- subset(x, grepl("\\+", strand))
x.tss <- subset(x.positive, grepl("tss", type))
x.utr <- subset(x.positive, grepl("5\\'utr", type))

x.sort <- x.tss[order(x.tss$pos),]

x.agg <-aggregate(x.sort, list(Region = x.sort$pos), mean)

x.narrow <- subset(x.sort, grepl("NarrowPeak", shape))
x.narrow_sub <- x.narrow[sample(nrow(x.narrow), 83464), ]
x.n.agg <-aggregate(x.narrow_sub, list(Region = x.narrow_sub$pos), mean)

x.broad <- subset(x.sort, grepl("BroadWithPeak", shape))
x.broad_sub <- x.broad[sample(nrow(x.broad), 83464), ]
x.b.agg <-aggregate(x.broad_sub, list(Region = x.broad_sub$pos), mean)

x.weak <- subset(x.sort, grepl("WeakPeak", shape))
x.weak_sub <- x.weak[sample(nrow(x.weak), 83464), ]
x.w.agg <-aggregate(x.weak_sub, list(Region = x.weak_sub$pos), mean)

#rw30 <- rep(1/101,101)
#x.sort$nuclAbunRW <- filter(x.sort$nuclAbun, rw30, sides=2)
#x.narrow$nuclAbunRW <- filter(x.narrow$nuclAbun, rw30, sides=2)
#x.broad$nuclAbunRW <- filter(x.broad$nuclAbun, rw30, sides=2)
#x.weak$nuclAbunRW <- filter(x.weak$nuclAbun, rw30, sides=2)

##ggplot2
stat_sum_single <- function(fun, geom="point") {
  stat_summary(fun.y=fun, colour="red", geom=geom, size = 0.1)
}

png("./all_ggplot2_agg_ribon7.png", width = 1920, height = 960, units = "px")
p <- qplot(pos, nuclAbun, data=x.narrow, ylim=c(0,25))

p + 
  stat_sum_single(mean, geom="line")#, colour="blue") +
  stat_sum_single(median, geom="line")#, colour="green") +
  stat_sum_single(mode, geom="line")#, colour="red")

#line plot overlays
p <- ggplot()
p + 
  geom_ribbon(data=x.agg, aes(x=pos, ymin=0, ymax=nuclAbun)) +
  geom_line(data=x.n.agg, aes(x=pos, y=nuclAbun), colour="blue") +
  geom_line(data=x.b.agg, aes(x=pos, y=nuclAbun), colour="red") +
  geom_line(data=x.w.agg, aes(x=pos, y=nuclAbun), colour="green")
dev.off()

##facet plots
png("./At_full_genome_normalised.png", width = 1920, height = 960, units = "px")
p <- ggplot() + theme(panel.grid.major = element_line(size = .5, color = "grey")) + coord_cartesian(ylim = c(0, 25), xlim = c(-500, 2000)) + geom_vline(xintercept = 0, colour="blue", size = 2, linetype="dotted")
p1 <- p + geom_ribbon(data=x.agg, aes(x=pos, ymin=0, ymax=nuclAbun)) + labs(title = "At genome - All Data", x="Distance from TSS modal peak", y="Nucleosome abundance count")
p2 <- p + geom_ribbon(data=x.n.agg, aes(x=pos, ymin=0, ymax=nuclAbun)) + labs(title = "At genome - Narrow TSS Sites", x="Distance from TSS modal peak", y="Nucleosome abundance count")
p3 <- p + geom_ribbon(data=x.b.agg, aes(x=pos, ymin=0, ymax=nuclAbun)) + labs(title = "At genome - Broad TSS Sites", x="Distance from TSS modal peak", y="Nucleosome abundance count")
p4 <- p + geom_ribbon(data=x.w.agg, aes(x=pos, ymin=0, ymax=nuclAbun)) + labs(title = "At genome - Weak TSS Sites", x="Distance from TSS modal peak", y="Nucleosome abundance count")

multiplot(p1,p2,p3,p4, cols=2)
dev.off()

