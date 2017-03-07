library("sleuth")

base_dir <- "/home/sbi6dap/Projects/TRI"

sample_id <- dir(file.path(base_dir, "kallisto"))
sample_id
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir,  "kallisto", id))
kal_dirs

s2c <- read.table(file.path(base_dir, "tri_map.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = sample, condition)
s2c

s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)

so <- sleuth_prep(s2c, ~ condition)
so<- sleuth_fit(so)
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
so <- sleuth_wt(so, 'conditionWT')
sleuth_live(so)


results_ordered <- results_table[order(results_table$qval),]
table(results_ordered$qval <= 0.05)
write.table( subset(results_ordered, qval <= 0.05), file='sleuth.DE_transcripts.qval_0.05.txt', sep="\t",row.names=F, quote=F)
