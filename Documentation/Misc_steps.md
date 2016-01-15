
<h3>Misceleneous stages</h3>

<b>Subsample BAM file</b>
```
cat <(samtools view -H ES12.bam) <(samtools view -q 255 ES12.bam | shuf -n 10000000) > ES12_10m.sam &
samtools view -bS ES12_10m.sam > ES12_10m.bam
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
