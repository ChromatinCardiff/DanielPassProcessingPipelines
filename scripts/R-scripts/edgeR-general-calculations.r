

x<- read.table("/home/sbi6dap/Projects/ALD/RNAseq/HTSeq/edgeR_analysis/toptags_edgeR_renamed.csv", header=TRUE, sep=",")
summary(x)

x$lightAvg <- ((x$ES1 + x$ES3 + x$ES5 + x$ES7)/4)
x$darkAvg <- ((x$ES2 + x$ES4 + x$ES6 + x$ES8)/4)

both.off <- subset(x, lightAvg == 0 & darkAvg == 0)
light.off.sig  <- subset(x, lightAvg == 0 & darkAvg != 0 & FDR <0.05)
dark.off.sig  <- subset(x, lightAvg != 0 & darkAvg == 0 & FDR <0.05)

light.upreg  <- subset(x, logFC <0 & FDR <0.05)    ## Negative numbers = light upregulated and vice versa
dark.upreg  <- subset(x, logFC >0 & FDR <0.05)

write.csv(light.upreg$Gene, file="light_upreg_sig_tags.csv")
write.csv(dark.upreg$Gene, file="dark_upreg_sig_tags.csv")
