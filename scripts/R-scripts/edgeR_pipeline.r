library("edgeR")
setwd("~/Projects/ALD/RNAseq/HTSeq")
samples <- read.table("samples.txt", header=TRUE)

counts = readDGE(samples$countf)$counts

noint = rownames(counts) %in% c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")
cpms = cpm(counts)
keep = rowSums(cpms >1) >=4 & !noint
counts = counts[keep,]

colnames(counts) = samples$shortname
head(counts[,order(samples$condition)], 5)

d = DGEList(counts=counts, group=samples$condition)
d = calcNormFactors(d)
plotMDS(d, labels=samples$shortname, col=c("darkgreen","blue", "darkred")[factor(samples$condition)])
legend("topleft", c("Dark (16h)", "Light (16h)", "Light (5d)"), fill=c("darkgreen","blue","darkred"))

d = estimateCommonDisp(d)
d = estimateTagwiseDisp(d)

plotMeanVar(d, show.tagwise.vars=TRUE, NBline=TRUE)
plotBCV(d)

de = exactTest(d, pair=c("16h-Light", "16h-Dark"))

tt = topTags(de, n=nrow(d))
head(tt$table)
nc = cpm(d, normalized.lib.sizes=TRUE)
rn = rownames(tt$table)
head(nc[rn,order(samples$condition)],5)

deg = rn[tt$table$FDR < .05]
plotSmear(d, de.tags=deg)
write.csv(tt$table, file="toptags_edgeR.csv")