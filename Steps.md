# DansProcessingPipeline

<h2>General steps that I've taken. Obviously varies case by case.</h1>

<h3>Bowtie command</h3>
```
bowtie -v 3 --trim3 14 --maxins 5000 --fr -k 1  -p 12 --sam indexes/TAIR10 -1 ES10_TCGGCA_L002_R1_001.fastq -2 ES10_TCGGCA_L002_R2_001.fastq bowtie1/ES10_50bp.sam -t
```

<h3>To parse SAM file into particle sizes</h3>
```
SAMparser2.pl -i infile.sam -o outdirectory -f HISEQ -p 0,80,150,350,500,680,860
```

<h3>To extract Chromosomes from sam file (not cleverly)</h3>
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
<h3>To generate histograms (for each file/chomosome)</h3>
```
for i in *_grep.txt; do (cut -f9 $i | sed 's/-//' | sort -n | uniq -c > $i.hist &); done
```
<h3>To merge chromosome files into one histogram (and addition of each chromosome values)</h3>
```
cat *.hist > temp.txt
awk '{arr[$2]+=$1} END {for (key in arr) printf("%s\t%s\n", key, arr[key])}' temp.txt   | sort +0n -1 >ES09_histogram.txt
```
Use basic_smooth-norm.r to chart

<h3>Job submission to slurm cluster</h3>
```
sbatch slurm_submission.sh
```
