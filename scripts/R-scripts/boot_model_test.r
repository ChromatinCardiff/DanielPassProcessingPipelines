library(boot)
library(plotrix)

DFlist <- list()

toload <- c(10, 14, 12)
for(i in toload)
{
  name <- paste("ES", i, ".xls", sep="")
  if(i == 9)
  {
    name <- paste("ES0", i, ".xls", sep="")
  }
  DFlist[[i]] <- read.table(paste("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/mean/mean_profile_TSS_heatmap/", name, sep=""), header=TRUE, sep="\t", row.names=1)
}

WT1 <- DFlist[[10]]
WT2 <- DFlist[[14]]
MU1 <- DFlist[[12]]
WT1B <- list()
WT2B <- list()

BootReps <- 200
StartMod <- -1500
EndMod <- 1490

meanBoot <- function(data, indices) {
	d <- data[indices]
	mean(d)
}

put.inNorm <- function(i, j, k) {
	pv <- pnorm(i, mean = j, sd = k)
	
	if(pv > 0.5)
	{
		pv <- 1 - pv
	}
	-log10(pv)
}

for(i in 1:ncol(WT1))
{
	if(i %% 5 == 0)
	{
		print(paste("Done: ", i, sep=" "))
	}
	WT1B[[i]] <- boot(WT1[,i], meanBoot, BootReps)
	WT2B[[i]] <- boot(WT2[,i], meanBoot, BootReps)
#	WT3B[[i]] <- boot(WT3[,i], meanBoot, BootReps)
}
allresults <- matrix(0, BootReps * 2, ncol(WT1))

for(i in 1:ncol(WT1))
{
	allresults[,i] <- c(WT1B[[i]]$t, WT2B[[i]]$t)
}

AllMeans <- apply(allresults, 2, function(x) mean(x))
AllSds <- apply(allresults, 2, function(x) sd(x))

muMeans <- colMeans(MU1)

# 'standardize' lines to zero mean and unit variance 
muMeans.d <- decostand(muMeans, 'standardize', MARGIN=2)
AllMeans.d <- decostand(AllMeans, 'standardize', MARGIN=2)

xNam <- seq(StartMod, EndMod, 50)

muPvs <- sapply(1:length(muMeans.d), function(i) put.inNorm(muMeans.d[i], AllMeans.d[i], AllSds[i]))

muPvs[is.infinite(muPvs)] <- 35

AllUpper <- AllMeans.d + (AllSds * 2)
AllLower <- AllMeans.d - (AllSds * 2)

print("Making graph...")

#pdf("BootsModel.pdf", width=15, height=10)
par(mfrow=c(2,1), mar=c(0,4,2,2))
par(fig=c(0,1,0.2,1))

plot(1:ncol(WT1), AllMeans.d, xlim=c(9,ncol(WT1)-8), type="l", ylab="Mean Read Counts", main="Distribution of Replicate Cumulative Counts vs Mutant", xaxt="n")
polygon(c(1:ncol(WT1), ncol(WT1):1), c(AllUpper, rev(AllLower)), density=NA, col="azure4")
points(1:ncol(WT1), AllMeans.d, type="l", lwd=1)
points(1:ncol(WT1), muMeans.d, type="l", lwd=3.5, col="firebrick1")
par(fig=c(0,1, 0, 0.2), new=TRUE, mar=c(3,4,0,2))

plot(1:ncol(WT1), muPvs, col="orange2", type="h", c(9, ncol(WT1)-8), lwd="3", ylab="Signif. -log(p)", xaxt = "n", xlab="Base Pairs Distance around TSS")
abline(h=-log10(5e-02))
abline(h=-log10((5e-02)/ncol(WT1)))
staxlab(1, at=seq(1, ncol(WT1), (ncol(WT1)/49.4)), labels=xNam, srt=45,  cex.axis=0.1, line.spacing=0.8, nlines=2, cex=0.7)

#dev.off()
