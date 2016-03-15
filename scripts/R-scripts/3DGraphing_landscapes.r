install.packages("plotly")
library(plotly)
library(MASS)
library(rgl)

options("scipen"=10000000000)

x1 <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/individual_genes/ES09-AT2G34420.3DG", header=TRUE)
x2 <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/individual_genes/ES10-AT2G34420.3DG", header=TRUE)
x3 <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/individual_genes/ES15-AT2G34420.3DG", header=TRUE)
x4 <- read.table("/home/sbi6dap/Projects/ALD/MNase-seq/individual_genes/ES16-AT2G34420.3DG", header=TRUE)

p <- ggplot(x, aes(alpha=0.2))
p + geom_point(aes(x=Pos, y=Sizes, colour=))

CDSlen <- seq(14522716,14523513)
genelen <- seq(14522524,14523568)


x1.den3d <- kde2d(x1$Pos, x1$Sizes, n= 500)
x2.den3d <- kde2d(x2$Pos, x2$Sizes, n= 500)
x3.den3d <- kde2d(x3$Pos, x3$Sizes, n= 500)
x4.den3d <- kde2d(x4$Pos, x4$Sizes, n= 500)

persp3d(x1.den3d$x, log(x1.den3d$y), x1.den3d$z, col=color[zcol], axes=TRUE, ticktype="detailed", zlab="")
lines3d(genelen,5,4e-06, col="black", lwd=25, lty=1)
lines3d(CDSlen,5,4e-06, col="blue", lwd=50, lty=1)
persp3d(x3.den3d$x, log(x3.den3d$y), x3.den3d$z, col=color[zcol], xlim=c(14522024,14524168), axes=TRUE, ticktype="detailed", zlab="")
persp3d(x3.den3d$x, log(x3.den3d$y), x3.den3d$z, col=color[zcol], xlim=c(14522024,14524168), axes=TRUE, ticktype="detailed", zlab="")
persp3d(x4.den3d$x, log(x4.den3d$y), x4.den3d$z, col=color[zcol], xlim=c(14522024,14524168), axes=TRUE, ticktype="detailed", smooth=FALSE)


persp3d(den3d$x, log(den3d$y), den3d$z, col=color[zcol], axes=TRUE, ticktype="detailed", zlab="")

nbcol = 100
color = rev(rainbow(nbcol, start = 0/6, end = 4/6))
zcol  = cut(den3d$z, nbcol)



plot_ly(x=den3d$x, y=den3d$y, z=den3d$z, type="surface") %>% layout(autosize = F, width = 1000, height = 100)

scene = list(xaxis=list(type = "log"))

plot_ly(den3d, x=x, y=y, z=z, type="surface") %>% layout(title="3d plot test", scene = scene)

