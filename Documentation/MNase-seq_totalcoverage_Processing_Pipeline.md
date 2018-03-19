<h3> Total Coverage processing</h3>
##Make sure you're using samtools version >1.2 otherwise things will break and you won't know why (well, now you know why)

<b>Map reads command</b>
```
bowtie -v 3 --trim3 14 --maxins 5000 --fr -k 1  -p 12 --sam indexes/TAIR10 -1 ES09_R1.fastq -2 ES09_R2.fastq bowtie1/ES09_50bp.sam -t
```
<b>convert mapped bam into paired-end bed file (coverage including gap between pairs), cut out relevant columns and sort</b>
```
samtools view -bf 0x2 ES09.bam | bamToBed -i stdin -bedpe | sort -k1,1 -k2,2n >ES09_bedpe.bed
#Or run with particle size filtering (here mapping only <120 length regions)
samtools view -bf 0x2 ES09.bam | bamToBed -i stdin -bedpe | awk '{if (($6 - $2) < 120) print $0}' | sort -k1,1 -k2,2n >ES09_lt120_bedpe.bed
```
<b>Calculate coverage of regions then contract into 10bp bins</b>
```
coverageBed -a chromosomes.bed -b ES09_bedpe.bed -d >ES09_totalcov.txt
cut -f1,4,5 ES09_totalcov.txt > ES09_totalcov_cut.txt
data_binner.pl -i ES09_totalcov_cut.txt -o ES09_totalcov.sgr -b 10 -m AVG
```

```
#Remove splice variant number from rownames
for i in *xls; do sed 's/^\(AT.......\)../\1/' $i >$i.cut; done
#Remove splice variant number from rownames AND duplicate row names
for i in *xls; do sed 's/^\(AT.......\)../\1/' $i | awk '!seen[$1]++' >$i.cut; done
```

# Run 2D vis for subNSP at TFBS

```
~/Projects/DansProcessingPipeline/scripts/core_scripts/BAM2SizePlot.py -i ES09.bam,ES13.bam -b ../REFDB/TFBS/PIF4.bed -S 100 -z None -q 1 -e 150
~/Projects/DansProcessingPipeline/scripts/core_scripts/BAM2SizePlot.py -i ES11.bam,ES15.bam -b ../REFDB/TFBS/PIF4.bed -S 100 -z None -q 1 -e 150
```
