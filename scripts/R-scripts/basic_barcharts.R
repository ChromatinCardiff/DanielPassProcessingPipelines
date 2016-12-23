library(ggplot2)
install.packages("reshape2")
library(reshape2)

setwd("/home/sbi6dap/Projects/ALD/totalcoverage/lt120/RAWDATA/results_individually_normalised/pooled/sums/")


x <- data.frame(read.table("allsums_gtupquart.txt", header=TRUE))
head(x)
y <- cast(x, smt_value ~ Sample)
vals <- x$smt_value
quantile(vals)

x <- read.csv("/home/sbi6dap/Projects/ACS/analysis/dpos/profile_TSS_heatmap/HO_GO_abundances.csv", header=TRUE)

x.melt <- melt(x, id=c("Keyword.Category","Functional.Category"))

x.split <- split(x.melt, x.melt$Keyword.Category)

summary(x.split[["GO Biological Process"]])

ggplot(x.split[["GO Molecular Function"]], aes(x=Functional.Category, y=value, fill=variable)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_brewer(palette="Dark2") +
  coord_flip()
  facet_wrap(~Keyword.Category, nrow=1)

p <- ggplot(data=x)
p + 
  geom_boxplot(aes(x=Sample, y=smt_value), outlier.shape=NA) +
  coord_cartesian(ylim = c(15,60))
