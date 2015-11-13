# Metagenomic assembly

# Align reads to reference
```
/opt/shared/group/smbpk/software/bowtie-1.1.1/bowtie -v 3 --maxins 1000 --fr -k 1 -p 64 --sam /home/GROUP-smbpk/sbi6dap/Projects/ALD/indexes/TAIR10 -1 /home/GROUP-smbpk/sbi6dap/Projects/Ath-reseq/Ath-reseq-R1.fastq -2 /home/GROUP-smbpk/sbi6dap/Projects/Ath-reseq/Ath-reseq-R2.fastq /home/GROUP-smbpk/sbi6dap/Projects/Ath-reseq/Ath-reseq.sam -t --un /home/GROUP-smbpk/sbi6dap/Projects/Ath-reseq/unmapped-Ath-reseq.sam
```
# Denovo assemble unmapped reads
```
/home/GROUP-smbpk/sbi6dap/localbin/megahit --kmin-1pass --k-min 27 --k-step 10 --k-max 127 -r unmapped-Ath-reseq_subset.fasta -o MH-unmapped-sub
```
# Annotate contigs greater than size(x)
```
/home/GROUP-smbpk/sbi6dap/localbin/blastn -db /home/GROUP-smbpk/sbi6dap/db/nt-nuc -query /home/GROUP-smbpk/sbi6dap/Projects/Ath-reseq/MH-unmapped-sub/final.contigs_1kb.fa -out /home/GROUP-smbpk/sbi6dap/Projects/Ath-reseq/MH-unmapped-sub/final.contigs_1kb_blast.txt -num_threads 32 -outfmt 7
```
# Extend contigs
```
/home/GROUP-smbpk/sbi6dap/localbin/SOURCE/SSPACE-STANDARD-3.0_linux-x86_64/SSPACE_Standard_v3.0.pl -l /home/GROUP-smbpk/sbi6dap/Projects/unmapped_extending/SSPACE_libfile.txt -s /home/GROUP-smbpk/sbi6dap/Projects/Ath-reseq/final.contigs.fa -x 0 -T 32 -b Unmapped_extending -p 1
```
# Close Gaps between contigs
```
GapCloser -b /home/GROUP-smbpk/sbi6dap/Projects/Ath-reseq/unmapped_extending/gapcloser_lib.txt -a /home/GROUP-smbpk/sbi6dap/Projects/Ath-reseq/unmapped_extending/unmapped_extending.final.scaffolds.fasta -o /home/GROUP-smbpk/sbi6dap/Projects/Ath-reseq/unmapped_extending/closed_gaps -t 32 -l 150
```
