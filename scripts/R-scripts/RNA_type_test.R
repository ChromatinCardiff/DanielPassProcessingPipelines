####################
## RNA type test  ##
####################
library(ggplot2)
library(reshape2)
gene.names <- read.table("geneIDs-name-type.txt", header=TRUE)   ### Format: ID \t gene_name \t RNA_type

full <- read.table("BMR_all.txt", header=TRUE)                   ### Format: ID \t Sample1 \t Sample2 \t Sample3
full.merge <- merge(full, gene.names, by="ID", all.x = TRUE)

full.melt <- melt(full.merge, na.rm=TRUE)

p <- ggplot(full.melt, aes(variable))
p + geom_bar(aes(weight = value, fill=type))
