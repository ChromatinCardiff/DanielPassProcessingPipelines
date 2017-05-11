# structuremax usage
install.packages("quantreg")
library("reshape2")

x <- as.matrix(read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/mean_ARA11_GMC_TSS_heatmap/ES14_filtered.xls", header=TRUE, check.names=FALSE, sep="\t", row.names=1))
HTseq <- read.table("/home/sbi6dap/Projects/ALD/RNAseq/ARA11/HTseq/alltags_edgeR.csv", header=TRUE, sep=",", row.names=1)

x_sub <- x[,100:200]

xmeans <- colMeans(x_sub[sample(nrow(x_sub),100),])
plot(xmeans, type="l")

find_peaks(xmeans)

find_peaks <- function (x, m = 15){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}
