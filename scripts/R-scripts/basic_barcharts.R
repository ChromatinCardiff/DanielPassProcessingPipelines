library(ggplot2)

x <- read.csv("/home/sbi6dap/Projects/ACS/analysis/dpos/profile_TSS_heatmap/HO_GO_abundances.csv", header=TRUE)

x.melt <- melt(x, id=c("Keyword.Category","Functional.Category"))

x.split <- split(x.melt, x.melt$Keyword.Category)

summary(x.split[["GO Biological Process"]])

ggplot(x.split[["GO Molecular Function"]], aes(x=Functional.Category, y=value, fill=variable)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_brewer(palette="Dark2") +
  coord_flip()
  facet_wrap(~Keyword.Category, nrow=1)
