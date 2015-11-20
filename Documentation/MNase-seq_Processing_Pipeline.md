
<h3>From FASTQ to SGR</h3>
<b>Bowtie command</b>
```
bowtie -v 3 --trim3 14 --maxins 5000 --fr -k 1  -p 12 --sam indexes/TAIR10 -1 ES10_TCGGCA_L002_R1_001.fastq -2 ES10_TCGGCA_L002_R2_001.fastq bowtie1/ES10_50bp.sam -t
```

<b>To parse SAM file into particle sizes</b>
```
SAMparser2.pl -i infile.sam -o outdirectory -f HISEQ -p 0,80,150,350,500,680,860
```

<b>To extract Chromosomes from sam file (not cleverly)</b>
```
chr_split.sh
	grep -w 'Chr1' ES16.sam >Chr1_grep.txt &
	grep -w 'Chr2' ES16.sam >Chr2_grep.txt &
	grep -w 'Chr3' ES16.sam >Chr3_grep.txt &
	grep -w 'Chr4' ES16.sam >Chr4_grep.txt &
	grep -w 'Chr5' ES16.sam >Chr5_grep.txt &
	grep -w 'mitochondria' ES16.sam >mito_grep.txt &
	grep -w 'chloroplast' ES16.sam >chloro_grep.txt &
```
<b>To generate histograms (for each file/chomosome)</b>
```
for i in *_grep.txt; do (cut -f9 $i | sed 's/-//' | sort -n | uniq -c > $i.hist &); done
```
<b>To merge chromosome files into one histogram (and addition of each chromosome values)</b>
```
cat *.hist > temp.txt
awk '{arr[$2]+=$1} END {for (key in arr) printf("%s\t%s\n", key, arr[key])}' temp.txt   | sort +0n -1 >ES09_histogram.txt
```
Use basic_smooth-norm.r to chart

<b>Generate SGR files</b>
```
sgr_builder.pl -i infile_150.txt -o outfile.sgr -p At_chr_sizes.txt

$ At_chr_sizes.txt
  Chr1	1251
	Chr2	55324
	Chr3	9876
	ABC		123
```

<h3>Using danpos</h3>

<b>convert sgr to wig</b>
```
sgr2wig.pl input.sgr output.wig
```
<b>Identify peaks</b>
```
#Numbers are reads which mapped to genome
danpos.py dpos ES09_150.wig,ES10_150.wig,ES11_150.wig,ES12_150.wig,ES13_150.wig,ES14_150.wig,ES15_150.wig,ES16_150.wig -o results_individually_normalised -c ES09_150.wig:112949625,ES10_150.wig:120774176,ES11_150.wig:107892210,ES12_150.wig:48206863,ES13_150.wig:75281419,ES14_150.wig:103625083,ES15_150.wig:95355763,ES16_150.wig:102368487
```
<b>Annotate positions by gene feature</b>
```
# Make sure chromosome labels EXACTLY match the genepred file! Case sensitive, fails with div0 error.
danpos.py profile results_individually_normalised/pooled/ES09_150.Fnor.smooth.wig,results_individually_normalised/pooled/ES10_150.Fnor.smooth.wig,results_individually_normalised/pooled/ES11_150.Fnor.smooth.wig,results_individually_normalised/pooled/ES12_150.Fnor.smooth.wig,results_individually_normalised/pooled/ES13_150.Fnor.smooth.wig,results_individually_normalised/pooled/ES14_150.Fnor.smooth.wig,results_individually_normalised/pooled/ES15_150.Fnor.smooth.wig,results_individually_normalised/pooled/ES16_150.Fnor.smooth.wig --genefile_paths ../../RNAseq-annotation-results/LDcommonmodel-ARA11.genepred --name mean_ARA11
```
<b>Annotate gtf file from cufflinks with genenames instead of generic labels and convert to genepred</b>
```
# To take out only the primary isoform
awk '{if ($6 == "1000") print $0;}' ES8-transcripts-expressed.gtf >ES8-primary-isoforms-only.gtf
# To convert to genepred
/home/sbi6dap/Projects/3rd_party_packages/UCSCtools/gtfToGenePred -genePredExt $i $i.genepred
rename 's/gtf.//' *
```
<b>danpos with RNA positions, only expressed genes</b>
danpos.py profile results_individually_normalised/pooled/ES11_150.Fnor.smooth.wig --genefile_paths ES8-primary-isoforms-only.genepred --flank_up 500

<b>join samples</b>
```
paste ES*TSS.xls |cut -f1,2,4,6,8,10,12,14,16 | sed 's/results_individually_normalised\/pooled\///g' | sed 's/.Fnor.smooth.wig.\/home\/sbi6dap\/Projects\/ALD\/RNAseq-annotation-results\//-/g' | sed 's/-primary-only.genepred//g' >all_TSS.csv
# for multiple RNA exp levels:
paste ES*TSS.xls |cut -f1,2,3,5,6,8,9,11,12 | sed 's/results_individually_normalised\/pooled\///g' | sed 's/.Fnor.smooth.wig.\/home\/sbi6dap\/Projects\/ALD\/RNAseq-annotation-results\/expression-split-test\/ES3_transcripts-FPKM-/-/g' | sed 's/.genepred.tss//g' >all_TSS.csv
```

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

#Remove splice variant number from rownames
```
for i in *xls; do sed 's/^\(AT.......\)../\1/' $i >$i.cut; done
```
