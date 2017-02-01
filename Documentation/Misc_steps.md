
<h3>Misceleneous stages</h3>

<b>Subsample BAM file</b>
```
cat <(samtools view -H ES12.bam) <(samtools view -q 255 ES12.bam | shuf -n 10000000) > ES12_10m.sam &
samtools view -bS ES12_10m.sam > ES12_10m.bam
samtools merge output.bam input1.bam input2.bam input*.bam
```
<b>Make insert size histogram</b>
```
cut -f9 AtN_rep.sam | grep -v '-' | sort | uniq -c >AtN_insertsize_hist.txt
```
<b>find closest gene feature</b>
```
#from bedops
closest-features --closest RAW_DATA_LightHigh-RAW_DATA_DarkHigh.positions.integrative.sort.bed ~/Projects/reference_data/AT_iGenome_genes.bed >test.bed
```
<b>Align on first nucleosome (or other feature)</b>
```
#Find max point in peak range (default: first 200bp after TSS) and slide trace to align.
plus1_slider.pl -i ES10_150.Fnor.smooth.wig.heatmap.xls -o ES10_150-1nucl-align

#Alternatively, output the shift values to a file for reuse later (i.e. calculate +1 nucleosome alignment for 150bp data and use the shifts for small particles)
plus1_slider.pl -i ES10_150.Fnor.smooth.wig.heatmap.xls -o ES10_150-1nucl-align -X 		#create file
plus1_slider.pl -i ES10_80.Fnor.smooth.wig.heatmap.xls -o ES10_80-1nucl-align -s ES10_150-1nucl-align-shift.txt  #use file instead of denovo calculations
```
<h3> Total Coverage processing</h3>
<b>convert total mapped bam into paired-end bed file (coverage including gap between pairs), cut out relevant columns and sort
```
samtools view -bf 0x2 ES09.bam | bamToBed -i stdin -bedpe | cut -f1,2,6 | sort -k1,1 -k2,2n >ES09_bedpe.txt
```
<b>Calculate coverage of regions then contract into 10bp bins</b>
```
coverageBed -a chromosomes.bed -b ES09_bedpe.bed -d >ES09_totalcov.txt
cut -f1,4,5 ES09_totalcov.txt > ES09_totalcov_cut.txt
data_binner.pl -i ES09_totalcov_cut.txt -o ES09_totalcov.sgr -b 10 -m AVG
```

```
#Remove splice variant number from rownames
for i in *DGE.xls; do sed 's/^\(AT.......\)\.[0-9]\+/\1/' $i >$i-gene; done
#Remove splice variant number from rownames AND duplicate row names
for i in *xls; do sed 's/^\(AT.......\)\.[0-9]\+/\1/' $i | awk '!seen[$1]++' >$i.cut; done
```

<b>Extract fasta sequence from a genome (bed format: chr \t start \t end). Can use gff instead of bed with same -bed flag.</b>
```
bedtools getfasta -fi Bos_taurus.UMD3.1.dna.toplevel.fa -bed tmp.bed -fo mir23a_27a_24-2_604up1757down.fasta
```

<b>Extract alignments that overlap with bed region</b>
```
intersectBed -abam BMR01out/accepted_hits.bam -b ../REFS/tmp.bed >ENSBTAG00000029996.bam
```

<b>extract gff columns to make bed format</b>
```
awk '{print $1 ":" $4 "-" $5 "\t" $0}' edit_Ler_cuffmatch_gene.gff >edit_Ler_cuffmatch_gene2.gff
```

<b>Coverage Histograms</b>(Ref: http://www.gettinggeneticsdone.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html)
```
bedtools coverage -hist -abam samp.01.bam -b target_regions.bed | grep ^all > samp.01.bam.hist.all.txt
```

<b>merge two files and fill gaps with zeros</b>
```
join -e '0' -o auto -a1 -a2 sort_uplist_GO.tsv sort_downlist_GO.tsv >updown_GO.tsv
```

<b>merge two files and numerically add (sum) a column</b>
```
paste -d+ <(cut -f5 prinegWT_tcv1.txt) <(cut -f5 prinegWT_tcv2.txt) | bc | paste <(cut -f1-4 prinegWT_tcv1.txt) - >tmp.txt
```

<b>Sort by first column only (-u uniq, -t delimiter, -k1,1 key field is col1)</b>
```
sort -u -t, -k1,1 file
```
<b>Pull out isoforms, reheader and drop danpos leading cols</b>
```
for i in *xls;
do
grep -f ~/Projects/AGM/RNAseq/cufflinks/isoforms/isoform_list.txt $i >tmp.txt;
cat HEADER.txt tmp.txt >tmp2.txt;
cut -f1,5- tmp2.txt >$i.cut;
done
```

<b>Subtract one wig from the other</b>
paste sample1.wig sample2.wig >tmp.wig
awk '{a=$1-$2; print a,"TABCHARACTER",$0}' tmp.wig | sed 's/^0\s\+fix/fix/' | cut -f1 >samples_sub.wig

<b>Using Homer to find TFBS locations</b>
```
# Create Motif file
seq2profile.pl <consensus> [# mismatches] [name] > output.motif
i.e. seq2profile.pl GGAAGT 0 ets > output.motif

# Find motif in genome
scanMotifGenomeWide.pl <file.motif> <genome.fasta> [-bed] > output.bed
scanMotifGenomeWide.pl TFBS/SORLREP2.motif TAIR10_Chr.all.fasta -bed >TFBS/SORLREP2.bed
```

# Convert file from Athamap into 10bp width bed:
awk '{$98 = $4-5; $99 = $4+5; print $1 ":" $98 "-" $99}' Athamap-PIF5.txt | sed 's/At\(.\).*:/Chr\1:/' > AthamapPIF5.bed
## Note: Clean up first and last lines!
