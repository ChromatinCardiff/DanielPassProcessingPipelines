library(vegan)
library(ggplot2)
library(ellipse)
library(grid)


# Load data file
x <-read.csv("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks/TSS_for_NMDS_500_upstream.txt", sep="\t", header=TRUE)
x <-read.csv("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks/TSS_for_NMDS_0-1000.txt", sep="\t", header=TRUE)
x.env <-read.csv("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks/sample_map.txt", sep="\t", header=TRUE)

#Transpose so its one rows:sample, column:OTU
rownames(x) <- x$pos
x = x[-1]
x_matrix <- t(x)

#Calculate distance
x.dis <- vegdist(x_matrix)

#View stresses
x.mds0 <-monoMDS(x.dis)
stressplot(x.mds0, x.dis)

#Make MDS
x.mds <- metaMDS(x_matrix, trace = FALSE)

plot(x.mds, type = "n")		##t = labeled, p = points
text(x.mds,display="sites",col="red",cex=1)
text(x.mds,display="species",col="black",cex=0.6)
fac1<-x.env$Treatment			## choose factor for grouping
fac2<-x.env$Exposure
fac3<-x.env$Pairs
ordispider(x.mds,group=fac1, label=TRUE, lwd=4, lty=3)
ordispider(x.mds,group=fac2, label=TRUE, lwd=4, lty=5)
ordispider(x.mds,group=fac3, label=FALSE, lwd=5, lty=1, col="blue")
ordiellipse(x.mds,group=fac1, label=TRUE, lwd=5)
ordiellipse(x.mds,group=fac2, label=TRUE, lwd=5)

#Extract NMDS coordinates
sites <- scores(x.mds, display = "sites")

#Copy coordinates out
rownames(x.env) <- x.env$Sample
x.env$Sample <-NULL
all <- merge(sites, x.env, by=0)

# Load in new table
#all <-read.csv("/home/sbi6dap/Projects/ALD/MNase-seq/dpos_peaks/sample_map.txt", sep="\t", header=TRUE)

df_ell <- data.frame()

for(g in levels(all$Exposure)){
df_ell <- rbind(df_ell, cbind(as.data.frame(with(all[all$Exposure==g,], ellipse(cor(NMDS1, NMDS2), 
        scale=c(sd(NMDS1),sd(NMDS2)), 
        centre=c(mean(NMDS1),mean(NMDS2))))),group=g))}

# Generate plot
 p <- ggplot(all, aes(NMDS1, NMDS2))

#plot it
 p + geom_point(size = 4, alpha=.8, aes(colour = factor(Treatment))) + geom_path(data=df_ell, aes(x=x, y=y, colour=group), size=1, linetype=1)
 p + geom_point(size = 4, alpha=.8, aes(colour = factor(Condition), shape = factor(Condition))) + geom_path(data=df_ell, aes(x=x, y=y, colour=group), size=1, linetype=1)


ef <- envfit(x.mds ~ Al_soil + P_soil + S_soil + Ca_soil + Ti_soil + V_soil + Cr_soil + Mn_soil + Fe_soil + Co_soil + Ni_soil + Cu_soil + Zn_soil + As_soil + Se_soil + Mo_soil + Cd_soil + Sb_soil + Ba_soil + Hg_soil + Pb_soil + pH_soil + moi_soil + OC_soil + SampNo_metal + Al_tissue + Mn_tissue + Fe_tissue + Co_tissue + Ni_tissue + Cu_tissue + Zn_tissue + Se_tissue + Sr_tissue + Mo_tissue + Cd_tissue + Sb_tissue + Pb_tissue + Hg_tissue + Cr_tissue + As_tissue, all)
ef.df<-as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r))		##
ef.df$species<-rownames(ef.df)						##


p + geom_point(size = 4, alpha=.8, aes(colour = factor(Species), shape = factor(Site))) + 
   geom_path(data=df_ell, aes(x=x, y=y, colour=group), size=1, linetype=1) + 
   geom_segment(data=ef.df,aes(x=0,xend=NMDS1,y=0,yend=NMDS2), arrow = arrow(length = unit(0.5, "cm")),colour="grey",inherit_aes=FALSE) + 
   geom_text(data=ef.df,aes(x=NMDS1,y=NMDS2,label=species),size=5)

## Examples: 	http://docs.ggplot2.org/current/geom_point.html

ord <-ordiellipse(x, x$worm, display ="sites", kind = "se", conf = 0.95, label = T)

## ellipses example: http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
## Ordinance types in R: http://www.pmc.ucsc.edu/~mclapham/Rtips/ordination.htm


x.sub0.150 <- x_matrix[c(0:150)]

PCoA.res<-capscale(x_matrix ~1,distance="euclidean")
plot(PCoA.res)
text(PCoA.res,display="sites",col="red",cex=0.6)
ordispider(PCoA.res,group=fac1, label=TRUE, lwd=1)
ordispider(PCoA.res,group=fac2, label=TRUE, lwd=1, lty=2)
ordiellipse(PCoA.res,group=fac2, label=TRUE, lwd=1, lty=2)
