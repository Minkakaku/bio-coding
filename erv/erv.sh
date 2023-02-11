# 下载参考组
find . -name '*' -exec basename {} .fq.gz \; > name

cat name | while read id;do
bsub -q normal_1week -n 12 -e $id.err -o $id.log -J $id \
trimgalore -w 12 -i 00.rawdata/$id/$id\_1.fq.gz -o QC/${id}.clean.fq.gz \
-j QC/$id.fastp.json -h QC/$id.fastp.html
done

cat name| while read id ; do
bsub -q normal_1week -n 12 -e $id.err -o $id.log -J $id \
hisat2 -t -x /data/gpfs02/wmo/genome/hisat2_index/hisat_offical_mm10/mm10/genome \
-U QC/${id}.clean.fq.gz -S 01.bam/$id.sam
done

cat name| while read id ; do
bsub -q normal_1week -n 20 -e $id.err -o $id.log -J $id \
samtools view -@ 20 -b 01.bam/$id.sam -o 01.bam/$id.bam
done

cat name|while read id;do
bsub -q normal_1week -n 8 -e $id.err -o $id.log -J $id \
samtools sort -@ 8 -o 01.bam/$id.sorted.bam 01.bam/$id.bam
samtools index -@ 8 01.bam/$id.sorted.bam
done

cat name|while read id;do
bsub -q normal_1week -n 8 -e $id.err -o $id.log -J $id \
samtools index -@ 8 01.bam/$id.sorted.bam
done

# counts
## 参考组有点问题，诺和的用v3，其他的用原版
cat name| while read id ; do
bsub -q normal_1week -n 12 -e $id.err -o $id.log -J $id \
featureCounts -T 12 -t exon -g gene_id -s 1 -a /data/gpfs02/wmo/genome/mm10/rep/mm10_rmsk_TE.gtf \
-o 02.count/$id.txt 01.bam/$id.sorted.bam
done

cat name| while read id ; do
bsub -q normal_1week -n 12 -e $id.err -o $id.log -J $id \
featureCounts -T 12 -t exon -g gene_id -s 1 -a /data/gpfs02/wmo/genome/mm10/rep/ucsc_mm10.rep_v3.gtf \
-o 02.count/$id.v3.count 01.bam/$id.sorted.bam
done

# cut
cat name| while read id;do
cut -f 1,7 02.count/$id.v3.count|grep -v '^#'> 02.count/$id.counts.txt
done

# legth
cat name| while read id;do
cut -f 1,6 02.count/$id.v3.count|grep -v '^#'> 02.count/$id.lengths.txt
done

#merge 
python /data/gpfs01/wmo/script/RNAseq/merge_metaphlan_tables.py 02.count/p60* > 02.count/p60_TE.counts.txt
