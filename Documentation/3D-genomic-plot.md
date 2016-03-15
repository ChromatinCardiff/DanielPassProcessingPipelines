
<h3>Steps to make continuous 3D genomic plot</h3>


```
# Co-ordinate bed file in format:
Chr2    14521524        14524568

# Extract bam region
intersectBed -abam ~/RawData/ALD/MNase-seq/SAMs/ES16.bam -b tmp.bed >ES16-AT2G34420.bam &
# Convert to sam
samtools view ES10-AT2G34420.bam >ES10-AT2G34420.sam
# Deconstruct co-ords and size
3DGenomic_size_plot.pl -i ES10-AT2G34420.sam -o ES10-AT2G34420.3DG

# input into 3DGenomic_landscape.r
```
