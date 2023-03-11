fasterq-dump -e 24 00.rawdata/$i/$i.sra -O 00.rawdata
trim_galore  -q 30 --phred33 -J 24 \
            --paired 00.rawdata/$i\_1.fastq 00.rawdata/$i\_2.fastq \
            --gzip -o 01.qc
hisat2 -t -p 24 -x /data/gpfs02/wmo/genome/hisat2_index/hisat2_mm39/mm39 \
-1 01.qc/$i\_1_val_1.fq.gz -2 01.qc/$i\_2_val_2.fq.gz -S 02.bam/$i.sam
samtools view -@ 24 -b 02.bam/$i.sam -o 02.bam/$i.bam
samtools sort -@ 24 -o 02.bam/$i.sorted.bam 02.bam/$i.bam
samtools index -@ 24 02.bam/$i.sorted.bam
featureCounts -T 24 -t exon -p -g gene_id -s 1 -a /data/gpfs02/wmo/genome/mm39/TE/mm39.TE.gtf \
-o 03.counts/$i.count 02.bam/$i.sorted.bam
cut -f 1,2,3,4,7 03.counts/$i.count |grep -v '^#' > 03.counts/$i.counts.txt
cut -f 1,6 03.count/$i.count |grep -v '^#'> 02.count/$id.lengths.txt