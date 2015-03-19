install.packages("ggplot2")
library("ggplot2")

x <- read.csv("/home/sbi6dap/Projects/ACS/analysis/dpos/peakcalling/dpeak/pooled/Athal_4_Part150_AB_10.smooth.positions.csv", header=TRUE)
colnames(x) <-c("chr", "start", "end", "smt_pos", "smt_value", "fuzziness_score")
#options(scipen=5)

x$logval <- log10(x$smt_value)

x.limit <- x[x$smt_value >=10 & x$smt_value <= 50,]

summary(x)

p <- ggplot(x.limit, aes(smt_value, fuzziness_score))
p + 
  geom_point(aes(colour = chr), alpha = 0.5, size = 1) +
  stat_smooth(method=lm, formula = y ~ x + I(x^2), aes(colour = factor(chr)))




png("./all_tss_line_rw100_inc_avg.png", width = 1920, height = 960, units = "px")


legend(1700,20, c("Narrow", "Broad", "Weak"), lty=c(1,1,1), lwd=c(2.5,2.5,2.5), cex=2, col=c("red", "green", "blue"))
dev.off()
