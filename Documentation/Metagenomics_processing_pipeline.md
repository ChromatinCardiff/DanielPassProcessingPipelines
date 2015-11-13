# Metagenomic assembly

<b>Align reads to reference</b>
```
/opt/shared/group/smbpk/software/bowtie-1.1.1/bowtie -v 3 --maxins 1000 --fr -k 1 -p 64 --sam /home/GROUP-smbpk/sbi6dap/Projects/ALD/indexes/TAIR10 -1 /home/GROUP-smbpk/sbi6dap/Projects/Ath-reseq/Ath-reseq-R1.fastq -2 /home/GROUP-smbpk/sbi6dap/Projects/Ath-reseq/Ath-reseq-R2.fastq /home/GROUP-smbpk/sbi6dap/Projects/Ath-reseq/Ath-reseq.sam -t --un /home/GROUP-smbpk/sbi6dap/Projects/Ath-reseq/unmapped-Ath-reseq.sam
```
<b>Denovo assemble unmapped reads</b>
```
/home/GROUP-smbpk/sbi6dap/localbin/megahit --kmin-1pass --k-min 27 --k-step 10 --k-max 127 -r unmapped-Ath-reseq_subset.fasta -o MH-unmapped-sub
```
<b>Extend contigs</b>
```
/home/GROUP-smbpk/sbi6dap/localbin/SOURCE/SSPACE-STANDARD-3.0_linux-x86_64/SSPACE_Standard_v3.0.pl -l /home/GROUP-smbpk/sbi6dap/Projects/unmapped_extending/SSPACE_libfile.txt -s /home/GROUP-smbpk/sbi6dap/Projects/Ath-reseq/final.contigs.fa -x 0 -T 32 -b Unmapped_extending -p 1
```
<b>Close Gaps between contigs</b>
```
GapCloser -b /home/GROUP-smbpk/sbi6dap/Projects/Ath-reseq/unmapped_extending/gapcloser_lib.txt -a /home/GROUP-smbpk/sbi6dap/Projects/Ath-reseq/unmapped_extending/unmapped_extending.final.scaffolds.fasta -o /home/GROUP-smbpk/sbi6dap/Projects/Ath-reseq/unmapped_extending/closed_gaps -t 32 -l 150
```

<b>Annotate contigs greater than size(x)</b>
```
~/shared_scripts/filter_short_fasta.pl 1000 <final.contigs.fa >final.contigs_1kb.fa

/home/GROUP-smbpk/sbi6dap/localbin/blastn -db /home/GROUP-smbpk/sbi6dap/db/nt-nuc -query /home/GROUP-smbpk/sbi6dap/Projects/Ath-reseq/MH-unmapped-sub/final.contigs_1kb.fa -out /home/GROUP-smbpk/sbi6dap/Projects/Ath-reseq/MH-unmapped-sub/final.contigs_1kb_blast.txt -num_threads 32 -outfmt 7
```

<b>Get taxonomic ID for matches</b>
```
grep 'hits found' -A1 final.contigs_1kb_blast.txt >final.contigs_1kb_blast_top.txt
grep -v '^#' final.contigs_1kb_blast_top.txt | cut -f1 | cut -f4 -d '|' >ACCs.txt

while read ACC ; do    echo -n -e "$ACC\t";    curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACC}&rettype=fasta&retmode=xml" |   grep TSeq_taxid |   cut -d '>' -f 2 |   cut -d '<' -f 1 |   tr -d "\n";    echo;  done <ACCs.txt >ACC_taxids.txt &
```
