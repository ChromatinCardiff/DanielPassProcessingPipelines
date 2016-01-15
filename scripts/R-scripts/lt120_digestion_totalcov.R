library(vegan)

setwd("/home/sbi6dap/Projects/ALD/totalcoverage/size_comparisons/ALL_TSS")
# Chromatin per-gene profile files
myFiles <- list.files(pattern="*DGE.xls-nodup")
filenames = c("ES09","ES11","ES13","ES15")

for(x in seq_along(myFiles)){
  data.in <- read.table(myFiles[x], header=TRUE, sep="\t", row.names=1)
  assign(paste(filenames[x]), data.in)
}