TSS150 <- read.table("/home/sbi6dap/Projects/AGM/MNaseseq/particles/150_norm/WholeGenome/WholeGenome_TSS.xls", header=TRUE, sep="\t", row.names=1)
TSS150_fil <- read.table("/home/sbi6dap/Projects/AGM/MNaseseq/particles/150_norm/WholeGenome-filter1_TSS.xls", header=TRUE, sep="\t", row.names=1)
TSSlt120 <- read.table("/home/sbi6dap/Projects/AGM/MNaseseq/lt120/lt120_norm/WholeGenome_ex1pc_TSS.xls", header=TRUE, sep="\t", row.names=1)
TSSgt120 <- read.table("/home/sbi6dap/Projects/AGM/MNaseseq/gt120/gt120_norm/WholeGenome_TSS.xls", header=TRUE, sep="\t", row.names=1)
TSSRS <- read.table("/home/sbi6dap/Projects/AGM/MNaseseq/RootVsShoot/150_norm/WholeGenome-ShootRef_TSS.xls", header=TRUE, sep="\t", row.names=1)
TSSRSgt120 <- read.table("/home/sbi6dap/Projects/AGM/MNaseseq/RootVsShoot/gt120_norm/WholeGenome-ShootRef_TSS.xls", header=TRUE, sep="\t", row.names=1)
TSSRSlt120 <- read.table("/home/sbi6dap/Projects/AGM/MNaseseq/RootVsShoot/lt120_norm/WholeGenome-ShootRef_TSS.xls", header=TRUE, sep="\t", row.names=1)


## quick tmp run
tmp <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks-RNA-guided/mean_ARA11_GMC_TSS.xls", header=TRUE, sep="\t", row.names=1)
tmpAnnot <- c("Light","Light","Dark","Dark","Light","Light", "Dark", "Dark")
tmpAnnot <- c("LightLow","LightHigh","DarkLow","DarkHigh","LightLow","LightHigh", "DarkLow", "DarkHigh")
plotLineGraph(tmp, tmpAnnot, "ALD TSS 150")


#List of annotations of the data
WOAnnot <- c("WT","WT","WT","OE","OE","OE")

# Give three factors: Data that you want to plot, the annotation list and a title for the graph
plotLineGraph(TSS150_ex, WOAnnot, "TSS - Whole Genome - 150bp")
plotLineGraph(TSSgt120, WOAnnot, "TSS - Whole Genome - gt120")
plotLineGraph(TSSlt120, WOAnnot, "TSS - Whole Genome - lt120")

#Annots
IAnnot <- c("Shoot1","Shoot3","Shoot5","Root1","Root2","Root3")
RSAnnot <- c("Shoot","Shoot","Shoot","Root-AGM34","Root","Root")

plotLineGraph(TSSRS, RSAnnot, "TSS - Whole Genome - 150bp")
plotLineGraph(TSSRSgt120, RSAnnot, "TSS - Whole Genome - gt120")
plotLineGraph(TSSRSlt120, RSAnnot, "TSS - Whole Genome - lt120")


plotLineGraph <- function(x, Annot, title){
  runData.dec <- as.data.frame(decostand(t(x), 'standardize', MARGIN=1))
  runData.dec$pos <- rownames(runData.dec)
  runData.dec = cbind(Annot,runData.dec)
  runData.melt <- melt(runData.dec, id=c("pos", "Annot"))

  # Chart all columns
  p <-ggplot(data=runData.melt, aes(x=as.numeric(as.character(variable)), y=value, colour=Annot))
  TSS.plot <-p +
    stat_smooth(method="loess", span=0.05, se=TRUE) +
    #geom_line(aes(x=as.numeric(as.character(variable)), y=value))+#, colour=Exposure)) +
    scale_x_continuous(breaks = pretty_breaks(n=12)) +
    geom_vline(xintercept=0, colour="blue", linetype="longdash") +
    theme(legend.position = "bottom") +
    ggtitle(title) +
    xlab("Base pairs from factor")
  
  return(TSS.plot)
  }

TSS150.plot
TSS150_ex.plot
TSS150_fil.plot
