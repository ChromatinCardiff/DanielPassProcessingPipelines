
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
# Numbers are reads which mapped to genome
danpos.py dpos ES09_150.wig,ES10_150.wig,ES11_150.wig,ES12_150.wig,ES13_150.wig,ES14_150.wig,ES15_150.wig,ES16_150.wig -o results_individually_normalised -c ES09_150.wig:112949625,ES10_150.wig:120774176,ES11_150.wig:107892210,ES12_150.wig:48206863,ES13_150.wig:75281419,ES14_150.wig:103625083,ES15_150.wig:95355763,ES16_150.wig:102368487
```
<b>Annotate positions by gene feature</b>
```
# Make sure chromosome labels EXACTLY match the genepred file! Case sensitive, fails with div0 error.
danpos.py profile results_individually_normalised/pooled/ES09_150.Fnor.smooth.wig,results_individually_normalised/pooled/ES10_150.Fnor.smooth.wig,results_individually_normalised/pooled/ES11_150.Fnor.smooth.wig,results_individually_normalised/pooled/ES12_150.Fnor.smooth.wig,results_individually_normalised/pooled/ES13_150.Fnor.smooth.wig,results_individually_normalised/pooled/ES14_150.Fnor.smooth.wig,results_individually_normalised/pooled/ES15_150.Fnor.smooth.wig,results_individually_normalised/pooled/ES16_150.Fnor.smooth.wig --genefile_paths ../../ACS/reference_data/at_tair10_mod.genepred --flank_up 500
```
