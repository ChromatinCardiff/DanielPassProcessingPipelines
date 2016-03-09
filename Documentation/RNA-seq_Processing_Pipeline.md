Some of the stages and commands I've been using. Not complete.
Lots of information can be found in http://www.nature.com/nprot/journal/v7/n3/full/nprot.2012.016.html

<h3>Mapping with Tophat</h3>
Get genome fasta and gff from https://www.ncbi.nlm.nih.gov/genome/ or similar, and turn gff into gtf with ```gffread file.gff -T -o file.gtf```. Or get gtf straight away, that'd work too.
```
/home/GROUP-smbpk/sbi6dap/localbin/tophat2 -p 32 -G /home/GROUP-smbpk/sbi6dap/working/RNA-seq/working/genes.gtf -o ES20_THout /home/GROUP-smbpk/sbi6dap/working/RNA-seq/working/genome fastqs/ES20_R1.fastq fastqs/ES20_R2.fastq
```

<h3>Transcriptome analysis cufflinks</h3>
```
/home/GROUP-smbpk/sbi6dap/localbin/cufflinks -g /home/GROUP-smbpk/sbi6dap/working/RNA-seq/working/genes.gtf -b /home/GROUP-smbpk/sbi6dap/working/RNA-seq/working/genome.fa -p 32 -o /home/GROUP-smbpk/sbi6dap/working/RNA-seq/working/CL-annot/ES8_CLout /home/GROUP-smbpk/sbi6dap/working/RNA-seq/working/TH/ES8_THout/accepted_hits_ES8.bam
```

<h4>Create reference gtf</h4>
```
cuffcompare -s genome.fa -CG -R -V -r genes.gtf CL-annot/ES5_CLout-Gtest/ES5_CLout-Gtest-transcripts.gtf CL-annot/ES2_CLout-Gtest/ES2_CLout-Gtest-transcripts.gtf CL-annot/ES6_CLout-Gtest/ES6_CLout-Gtest-transcripts.gtf CL-annot/ES3_CLout-Gtest/ES3_CLout-Gtest-transcripts.gtf CL-annot/ES7_CLout-Gtest/ES7_CLout-Gtest-transcripts.gtf CL-annot/ES1_CLout-Gtest/ES1_CLout-Gtest-transcripts.gtf CL-annot/ES4_CLout-Gtest/ES4_CLout-Gtest-transcripts.gtf CL-annot/ES8_CLout-Gtest/ES8_CLout-Gtest-transcripts.gtf

# Put gene_name as gene_id instead of "XLOC_1234567"

```

<h3>TSS prediction (aka pulling mapped primary isoform start sites)</h3>
```
awk '{if ($3 ~ /transcript/) {print $0;} }' transcripts.gtf >allgenes.gtf         #Pull out transcripts only, ignore individual exons
grep -v 'FPKM "0.0000000000"' allgenes.gtf >isoforms.gtf                          #Remove unexpressed lines
awk '{if ($6 ~ /1000/) {print $0;} }' isoforms.gtf >primary-isoforms-only.gtf       #extract only primary isoform (having both would be useful in some cases)
```

<h3>Differential expression analysis with HTseq and EdgeR</h3>

HTSeq, repeat for each sample
```
htseq-count -f 'bam' -s no -a 10 /home/GROUP-smbpk/sbi6dap/working/RNA-seq/working/ARA11/tophat/ES1_THout/accepted_hits.bam /home/GROUP-smbpk/sbi6dap/working/RNA-seq/Ath_ref_files/ARA11/Araport11_genes.20150701.gtf > /home/GROUP-smbpk/sbi6dap/working/RNA-seq/working/ARA11/HTSeq_out/ES1.count
```
EdgeR
```
library("edgeR")
setwd("~/Projects/ALD/RNAseq/HTSeq")
samples <- read.table("samples.txt", header=TRUE)

counts = readDGE(samples$countf)$counts

noint = rownames(counts) %in% c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")
cpms = cpm(counts)
keep = rowSums(cpms >1) >=4 & !noint
counts = counts[keep,]

colnames(counts) = samples$shortname
head(counts[,order(samples$condition)], 5)

d = DGEList(counts=counts, group=samples$condition)
d = calcNormFactors(d)
plotMDS(d, labels=samples$shortname, col=c("darkgreen","blue", "darkred")[factor(samples$condition)])
legend("topleft", c("Dark (16h)", "Light (16h)", "Light (5d)"), fill=c("darkgreen","blue","darkred"))

d = estimateCommonDisp(d)
d = estimateTagwiseDisp(d)

plotMeanVar(d, show.tagwise.vars=TRUE, NBline=TRUE)
plotBCV(d)

de = exactTest(d, pair=c("16h-Light", "16h-Dark"))

tt = topTags(de, n=nrow(d))
head(tt$table)
nc = cpm(d, normalized.lib.sizes=TRUE)
rn = rownames(tt$table)
head(nc[rn,order(samples$condition)],5)

deg = rn[tt$table$FDR < .05]
plotSmear(d, de.tags=deg)
write.csv(tt$table, file="toptags_edgeR.csv")
```
