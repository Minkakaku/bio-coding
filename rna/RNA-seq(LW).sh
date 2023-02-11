## LW RNA-seq(mouse)

## translate data to hpc
scp -r p53d* wmo@hpc.xmu.edu.cn:/data/gpfs02/wmo/work/HF/jobs/ZHRGSE150714/00.rawdata

cd /data/gpfs02/wmo/work/RNA-seq/LW_P53/
## QC
cat name|while read id
do
bsub -q normal -n 10 -e $id.err -o $id.log -J $id \
fastqc -t 10 01.cleandata/$id/$id\_1.fq.gz -o 00.QC

## Align
cat name|while read id;do
bsub -q normal_1day -n 24 -e $id.err -o $id.log -J $id \
STAR --runThreadN 24 --readFilesCommand zcat --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate \
--genomeDir /data/gpfs02/wmo/genome/reference/STAR_index_mm10_50bp --readFilesIn 01.cleandata/$id/$id\_1.fq.gz \
--outFileNamePrefix 02.Alignment/$id #^N

cat name|while read id;do
bsub -q normal_1week -n 8 -e $id.err -o $id.log -J $id \
samtools sort -@ 8 -o 02.Alignment/$id.sorted.bam 02.Alignment/$id\Aligned.sortedByCoord.out.bam #^N
samtools index -@ 8 02.Alignment/$id.sorted.bam #^N


## counting 
cat name|while read id;do
bsub -q normal -n 20 -e $id.err -o $id.log -J $id \
featureCounts -T 20 -p -t exon -g gene_id -a /data/gpfs01/wmo/genome/mm10/mm10.ucsc.gtf -o 04.featureCounts/$id.featureCounts.txt 03.Alignment/$id.sorted.bam
done
cat name|while read id;do

cut -f 1,7 04.Counts/$id.Counts.txt|grep -v '^#'>04.Counts/$id\.Counts.filter.tsv
done
#merge 
python ~/script/RNAseq/merge_metaphlan_tables.py 04.featureCounts/*[0-9].txt > 04.featureCounts/merge.all_RNA_counts.txt


