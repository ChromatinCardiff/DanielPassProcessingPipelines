x <- read.table("~/working/ACS/ACS/Original_data/2012_2013Prelim_Ang_Marg/Exp1/split_chr/chr4_3peat.txt", header=FALSE)
colnames(x) <-c("chr", "nuclAbun", "pos", "strand", "rnaCount", "shape", "type", "gene")
#options(scipen=5)

#x.not0 <- subset(x, !grepl(0, val))

x.positive <- subset(x, grepl("\\+", strand))
x.tss <- subset(x.positive, grepl("tss", type))

x.sort <- x.tss[order(x.tss$pos),]

x.narrow <- subset(x.sort, grepl("NarrowPeak", shape))
x.broad <- subset(x.sort, grepl("BroadWithPeak", shape))
x.weak <- subset(x.sort, grepl("WeakPeak", shape))

rw30 <- rep(1/101,101)
x.sort$nuclAbunRW <- filter(x.sort$nuclAbun, rw30, sides=2)
x.narrow$nuclAbunRW <- filter(x.narrow$nuclAbun, rw30, sides=2)
x.broad$nuclAbunRW <- filter(x.broad$nuclAbun, rw30, sides=2)
x.weak$nuclAbunRW <- filter(x.weak$nuclAbun, rw30, sides=2)

png("./all_tss_line_rw100_inc_avg.png", width = 1920, height = 960, units = "px")
plot(nuclAbun ~ pos, data=x.sort, ylim=c(0,20), type="n")
lines(x.sort$pos, x.sort$nuclAbunRW,  col="black", lwd=5, lty=2)
lines(x.weak$pos, x.weak$nuclAbunRW,  col="blue")
lines(x.broad$pos, x.broad$nuclAbunRW,  col="green")
lines(x.narrow$pos, x.narrow$nuclAbunRW,  col="red")
legend(1700,20, c("Narrow", "Broad", "Weak"), lty=c(1,1,1), lwd=c(2.5,2.5,2.5), cex=2, col=c("red", "green", "blue"))
dev.off()
