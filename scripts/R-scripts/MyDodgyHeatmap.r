library(ggplot2)
library(reshape)
library(dplyr)
library(RColorBrewer)
library(gplots)

#library(gplots)
IDtoGene <- read.delim("/home/daniel/GRH/YSR/Homo_sapiens.GRCh38.97.map.txt", header=T)


# Pick input
x = read.csv("/home/daniel/Dropbox/Work/GRH/Yasir_Syed/gProfiler/gProfiler_rmout_summary.csv", header = TRUE)
x = read.csv("/home/daniel/GRH/YSR/rmdup/ASD_norm_counts_sigonly.txt", header = TRUE)
x = read.csv("/home/daniel/GRH/YSR/rmdup/SCZ_norm_counts.txt", header = TRUE)
x = read.csv("/home/daniel/GRH/YSR/rmdup/Anysigs.csv", header = TRUE)


dat.annot <- merge(x, IDtoGene, by.x="Id", by.y="Id" )

rownames(dat.annot) <- dat.annot$Name 
dat.annot$Id <- NULL
dat.annot$Name <- NULL

#filt <- "TF"
#x.filt <- subset(x, source == filt)
#x.filt.cols <- select(x.filt, term_name, padj.DelUp, padj.DupUp, padj.DelDown, padj.DupDown)
#rownames(x.filt.cols) <- x.filt.cols$term_name
#x.filt.heatmap <- data.matrix(x.filt.cols[2:5])

my_palette <- colorRampPalette(brewer.pal(9,"Reds"))

heatmap.2(data.matrix(log2(dat.annot+1)),
          Rowv=TRUE,
          Colv=FALSE,                              #       if you enable dendrograms then these have to true
          #distfun = dist,
          #scale="column",
          #margins = c(7, 15),
          col = my_palette,
          cexRow = 0.65,
          cexCol = 0.95,
          #labRow = NULL,
          #labCol = NULL,
          sepcolor = "black",
          sepwidth=c(0.001,0.001),
#          colsep = 1:ncol(x.filt),
          #rowsep = 1:nrow(x.filt.heatmap),
          xlab = "samples", ylab = "genes", main = "SCZSigGenes",
          #key=TRUE,
          #keysize=1,
          trace="none",
          density.info=c("none")
)

x.melt <- melt(x.filt)

base_size <- 9
p <- ggplot(x.melt, aes(x = X2, y = X1)) + 
  geom_tile(aes(fill = log(value+1)), colour = "white") + 
  scale_fill_gradient(low = "steelblue", high = "white")

p + theme_grey(base_size = base_size) + labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0)) +  scale_y_discrete(expand = c(0, 0))


