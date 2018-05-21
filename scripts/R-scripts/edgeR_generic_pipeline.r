library(edgeR)
setwd("~/Projects/ALD/RNAseq/ARA11/HTSeq/")

# [optional] Pre-calculate length with gtf2lengthGC.r
LengthGC <- read.table("~/Projects/REFDB/GC_lengths.tsv", header=TRUE)

samples <- read.table("samples_cycd.txt", header=TRUE)

counts = data.frame(readDGE(samples$countfile)$counts)
noint = rownames(counts) %in% c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")
cpms = cpm(counts)
# [optional] requires the length calculations
rpkm <- rpkm(counts, LengthGC$Length)     # Warning message may appear because of some mismatched genes. Dont worry
rpkm.colsums <- colSums(rpkm)
TPM <- sweep(rpkm, 2, rpkm.colsums, `/`)
TPM <- apply(TPM, 2, function(x) (x * 10^6 ))

# Filter out <1 counts per million
keep = rowSums(TPM >3) >=3 & !noint  # 3 for smallest rep group
exp.matrix = TPM[keep,]
#ALTERNATIVE#
keep = rowSums(rpkm >3) >=3 & !noint  # 3 for smallest rep group
exp.matrix = rpkm[keep,]
#ALTERNATIVE#
keep = rowSums(cpms >3) >=3 & !noint  # 3 for smallest rep group
exp.matrix = counts[keep,]

## OR dont filter anything ##
nrow(exp.matrix)

colnames(exp.matrix) = samples$SampleID
head(exp.matrix[,order(samples$Modification)], 5)

d = DGEList(counts=exp.matrix, group=samples$Modification)        # Test by condition
d = calcNormFactors(d)

par(mar=c(5.1, 4.1, 4.1, 9.1), xpd=TRUE)
plotMDS(d, labels=samples$shortname, col=c("red", "blue")[factor(samples$Modification)])
legend("right", inset=c(-0.25,0), c("CYCD31-OE", "WT"), fill=c("red", "blue"))
title(main="MDS plot of sample relationships (paired)")
par(xpd=FALSE)

d = estimateCommonDisp(d)
d = estimateTagwiseDisp(d)

plotMeanVar(d, show.tagwise.vars=TRUE, NBline=TRUE)
title(main="Per-gene Mean-Variance relationship")
plotBCV(d)
title(main="biological coefficient of variation versus counts")

###############
# Comparisons #
###############
gene.names <- read.table("~/Projects/REFDB/geneIDs-name.txt", header=TRUE)

de1 = exactTest(d, pair=c("WT", "cycd"))        ## AKA Fischers exact test

tt = topTags(de1, n=nrow(d))
nc = cpm(d, normalized.lib.sizes=TRUE)
rn = rownames(tt$table)
deg = rn[tt$table$FDR < .05]
plotSmear(d, de.tags=deg)
title(main="WT vs CYCD31-OE")
tt$table$ID <- row.names(tt$table)
tmp.merge <- merge(tt$table, gene.names, by="ID", all.x = TRUE)
de1.sort <- tmp.merge[with(tmp.merge, order(FDR)),]
#de1.sort <- tt$table[with(tt$table, order(FDR)),]
head(de1.sort, n=10)
write.csv(de1.sort, file="toptags_edgeR-WT_vs_cycd.csv", row.names = FALSE)
