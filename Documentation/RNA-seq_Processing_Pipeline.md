Some of the stages and commands I've been using. Not complete.
Lots of information can be found in http://www.nature.com/nprot/journal/v7/n3/full/nprot.2012.016.html

<h3>Mapping with Tophat</h3>
Get genome fasta and gff from https://www.ncbi.nlm.nih.gov/genome/ or similar, and turn gff into gtf with `gffread file.gff -T -o file.gtf`. Or get gtf straight away, that'd work too.

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

<h3>TSS isoform selection (aka pulling mapped primary isoform start sites)</h3>

```
#Pull out transcripts only, ignore individual exons | Remove unexpressed lines | extract only primary isoform (having both would be useful in some cases)
for i in AGM*/transcripts.gtf;
do
awk '{if ($3 ~/transcript/) {print $0;} }' $i | grep -v 'FPKM "0.0' | awk '{if ($6 ~ /1000/) {print $0;} }' >$i-primary.gtf;
done

# Grab basic list of genes
for i in AGM*/*primary.gtf;
do
cut -f2 -d ';' AGM16/transcripts.gtf-primary.gtf | sed 's/ transcript_id "\(.*\)"/\1/' | grep 'AT' >$i-isoform_list.txt;
done

# Identify common isoforms between comparative samples
for i in AG*; do rename "s/tran/$i-tran/" $i/tra*; done
```

<h2>Annotating novel genome with RNAseq data</h2>

```
#Map to genome without reference gtf
tophat2 -p 32 -o AGM1_THout REFS/Ler0_FJ fastqs/AGM1_R1.fastq fastqs/AGM1_R2.fastq

# Assemble mapped reads
cufflinks -p 16 -o cufflinks/AGM19 AGM19_acc.bam

# Merge replicates
cuffmerge -s REFS/Ler0_FJ.fa -o WT_cuffmerge -p 8 WT_assembly_list.txt
cuffmerge -s REFS/Ler0_FJ.fa -o g_cuffmerge -p 8 g_assembly_list.txt

#Compare between treatments and make one final GTF
cuffcompare -o WT_vs_g_cuffcompare -s REFS/Ler0_FJ.fa WT_cuffmerge/merged.gtf g_cuffmerge/merged.gtf

# Extract fasta sequence for annotation
gffread WT_vs_g54_cuffcompare.combined.gtf -g ~/axon_scratch/REFS/Ler0_FJ.fa -w WT_vs_g54_cuffcompare.combined.fasta -i 450

# Annotate fasta
blastn -query WT_vs_g54_cuffcompare.combined.fasta -db ~/db/nt-nuc -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 16 -out WT_vs_g54_cuffcompare.combined.blastn
```

<h3>Make Coverage Chart for Genome Browser</h3>
```
for i in *bam; do bamToBed -i $i > $i.bed; done
for i in *bam.bed; do genomeCoverageBed -i $i -bg -g ~/Projects/DansProcessingPipeline/scripts/core_scripts/Atha_chr_sizes.txt > $i.cov; done
for i in *cov; do bedGraphToBigWig $i ~/Projects/DansProcessingPipeline/scripts/core_scripts/Atha_chr_sizes.txt $i.bigwig; done

# Or cat and sort Bed files together to make one representative trace, then make bigwig
cat ES1.bam.bed ES3.bam.bed ES5.bam.bed ES7.bam.bed | sort -k1,1 -k2,2 >ES-All_light.bed
```

<h3>Differential expression analysis with HTseq and EdgeR</h3>

HTSeq, repeat for each sample

```
htseq-count -f 'bam' -s no -a 10 /home/GROUP-smbpk/sbi6dap/working/RNA-seq/working/ARA11/tophat/ES1_THout/accepted_hits.bam /home/GROUP-smbpk/sbi6dap/working/RNA-seq/Ath_ref_files/ARA11/Araport11_genes.20150701.gtf > /home/GROUP-smbpk/sbi6dap/working/RNA-seq/working/ARA11/HTSeq_out/ES1.count

# Make single table in bash rather than touching R
cut -f1 tri1.count > all_tri.count
for i in tri*count
 do
  paste all_tri.count <(cut -f2 $i) > tmp.txt
  mv tmp.txt all_tri.count
 done
```
EdgeR

```
#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")
library(edgeR)
setwd("~/Projects/NON-ATHAL/BMR/HTSeq/")
# Pre-calculate length with gtf2lengthGC.r
LengthGC <- read.table("../REFS/GC_lengths.tsv", header=TRUE)

samples <- read.table("tmp_samples.txt", header=TRUE)

counts = data.frame(readDGE(samples$countfile)$counts)
noint = rownames(counts) %in% c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")
cpms = cpm(counts)
rpkm <- rpkm(counts, LengthGC$Length)
rpkm.colsums <- colSums(rpkm)
TPM <- sweep(rpkm, 2, rpkm.colsums, `/`)
TPM <- apply(TPM, 2, function(x) (x * 10^6 ))

# Filter out <1 counts per million                                  ############################
keep = rowSums(TPM >10) >=2 & !noint  # 2 for smallest rep group    #### Hash to not filter ####
exp.matrix = TPM[keep,]                                             ############################
#ALTERNATIVE#
keep = rowSums(rpkm >1) >=2 & !noint  # 2 for smallest rep group
exp.matrix = rpkm[keep,]
#ALTERNATIVE#
keep = rowSums(cpms >1) >=2 & !noint  # 2 for smallest rep group
exp.matrix = counts[keep,]
# OR dont filter anything
nrow(exp.matrix)

colnames(exp.matrix) = samples$origID
head(exp.matrix[,order(samples$condition)], 5)

d = DGEList(counts=exp.matrix, group=samples$condition2)        # Test by condition
d = calcNormFactors(d)

par(mar=c(5.1, 4.1, 4.1, 9.1), xpd=TRUE)
plotMDS(d, labels=samples$shortname, col=c("yellow", "red", "blue", "green", "darkgreen")[factor(samples$condition)])
legend("right", inset=c(-0.25,0), c("Ctrl-untreated", "Ctrl-Scrambled", "miR5p", "221", "222"), fill=c("darkgreen", "green", "yellow", "red", "blue"))
title(main="MDS plot of sample relationships (paired)")
par(xpd=FALSE)

d = estimateCommonDisp(d)
d = estimateTagwiseDisp(d)

plotMeanVar(d, show.tagwise.vars=TRUE, NBline=TRUE)
title(main="Per-gene Mean-Variance relationship")
plotBCV(d)
title(main="biological coefficient of variation versus TPM")

# Test by repeat number (for QC because of group 3 batching)
e = DGEList(counts=counts, group=samples$"repeat")
e = calcNormFactors(e)
e = estimateCommonDisp(e)
e = estimateTagwiseDisp(e)

###############
# Comparisons #
###############
gene.names <- read.table("geneIDs-name-type.txt", header=TRUE)

de1 = exactTest(d, pair=c("Untreated", "Neg-ctrl"))        ## AKA Fischers exact test
de2 = exactTest(d, pair=c("Neg-ctrl", "miR5p"))        ## AKA Fischers exact test

setwd("~/temp/analysis/")
pdf( "smearplots.pdf" , width = 10 , height = 7 )

tt = topTags(de1a, n=nrow(d))
nc = cpm(d, normalized.lib.sizes=TRUE)
rn = rownames(tt$table)
deg = rn[tt$table$FDR < .05]
plotSmear(d, de.tags=deg)
title(main="untreated_vs_neg")
tt$table$ID <- row.names(tt$table)
tmp.merge <- merge(tt$table, gene.names, by="ID", all.x = TRUE)
de1.sort <- tmp.merge[with(tmp.merge, order(FDR)),]
head(de1.sort)
write.csv(de1.sort, file="toptags_edgeR-untreated_vs_neg.csv", row.names = FALSE)
dev.off()
```
