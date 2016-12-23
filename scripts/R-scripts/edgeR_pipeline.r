library("edgeR")
setwd("~/Projects/AGM/RNAseq/HTSeq")

samples <- read.table("samples_LD.txt", header=TRUE)

counts = readDGE(samples$countf)$counts

noint = rownames(counts) %in% c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")
cpms = cpm(counts)

# Filter out <1 counts per million          ##############################
keep = rowSums(cpms >1) >=4 & !noint        #### ONLY DO IF FILTERING ####
counts = counts[keep,]                      #### ONLY DO IF FILTERING ####
# OR dont filter anything                   ##############################

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

de = exactTest(d, pair=c("16h-Light", "16h-Dark"))        ## AKA Fischers exact test

tt = topTags(de, n=nrow(d))
head(tt$table)
nc = cpm(d, normalized.lib.sizes=TRUE)
rn = rownames(tt$table)
head(nc[rn,order(samples$condition)],5)

deg = rn[tt$table$FDR < .05]
plotSmear(d, de.tags=deg)
write.csv(tt$table, file="alltags_edgeR-filter<1m.csv")


####### Gene barchart

library(ggplot2)
library(reshape2)

x <- data.frame(exp=d$pseudo.counts["AT4G12800", ])
#x <- data.frame(head(d$pseudo.counts, 10))
cats <- c("Light", "Dark", "Light", "Dark", "Light", "Dark", "Light", "Dark")
x <- cbind(x, cats)
x2 <- summarySE(x, measurevar="exp", groupvars=c("cats"))

ggplot(x2, aes(x=cats, y=log(exp))) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=log(exp-ci), ymax=log(exp+ci)), width=.2, position=position_dodge(.9))
  ylim(6.25,6.75)
  
facet_wrap(~Keyword.Category, nrow=1)

