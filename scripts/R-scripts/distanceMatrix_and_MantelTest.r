install.packages("ade4")
library(ade4)


sites <- read.table("/media/sbi6dap/HDDHome/Dropbox/Work/Projects/LRP/metadata/Site_coordinates.csv", header=TRUE,  check.names = FALSE)
unifrac <-  read.table("/media/sbi6dap/HDDHome/Dropbox/Work/Projects/LRP/Analysis/analysis/beta_diversity/weighted_unifrac_dm.txt", header=TRUE,  check.names = FALSE)

row.names(sites) <- sites$SampleID
sites.dm <- dist(cbind(sites$Lat, sites$Long))
as.matrix(sites.dm)[1:20, 1:20]

unifrac.dm <- as.dist(unifrac)
as.matrix(unifrac.dm)[1:20, 1:20]


mantel.rtest(sites.dm, unifrac.dm, nrepet = 9999)
