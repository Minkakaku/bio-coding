
cd /data/gpfs02/wmo/work/RNA-seq/ZBX_E18

mkdir 00.rawdata/ 01.fastq/ 02.cleandata/ 03.Align/
## download               
bsub -q normal_1week -n 1 -e down.err -o down.log -J down \
prefetch --option-file name -O 00.rawdata/


## translate srr file to fastq file
cat name |while read id;do
bsub -q normal -n 1 -e $id.err -o $id.log -J $id \
fastq-dump 00.rawdata/$id/$id.sra --gzip --split-e --defline-seq '@$ac-$si/$ri' --defline-qual '+' \
-O 01.fastq/
done

## QC
head name -n 14|while read id;do
bsub -q normal_1day_new -n 10 -e $id.err -o $id.log -J $id \
fastqc -t 10 01.fastq/$id.fastq.gz

tail name -n 14|while read id;do
bsub -q normal_1day_new -n 10 -e $id.err -o $id.log -J $id \
fastqc -t 10 01.fastq/$id\_1.fastq.gz
fastqc -t 10 01.fastq/$id\_2.fastq.gz


## trim
head name -n 14|while read id;do
bsub -q normal_1day_new -n 10 -e $id.err -o $id.log -J $id \
trim_galore -j 10 -q 20 --phred33 --stringency 3 --length 20 --fastqc -e 0.1 01.fastq/${id}.fastq.gz --gzip \
-o 02.cleandata/

tail name -n 14|while read id;do
bsub -q normal_1day_new -n 10 -e $id.err -o $id.log -J $id \
trim_galore -j 10 -q 20 --phred33 --stringency 3 --length 20 --fastqc --paired -e 0.1 01.fastq/${id}\_1.fastq.gz 01.fastq/${id}\_2.fastq.gz \
--gzip -o 02.cleandata/


## Align
head name -n 14|while read id;do
bsub -q normal_1week -n 24 -e $id.err -o $id.log -J $id \
STAR --runThreadN 24 --readFilesCommand zcat --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate \
--genomeDir /data/gpfs02/wmo/genome/reference/STAR_index_mm10_50bp --readFilesIn 02.cleandata/$id\_trimmed.fq.gz \
--outFileNamePrefix 03.Align/$id #^N

tail name -n 14|while read id;do
bsub -q normal_1week -n 24 -e $id.err -o $id.log -J $id \
STAR --runThreadN 24 --readFilesCommand zcat --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate \
--genomeDir /data/gpfs02/wmo/genome/reference/STAR_index_mm10_50bp --readFilesIn 02.cleandata/$id\_1_val_1.fq.gz 02.cleandata/$id\_2_val_2.fq.gz \
--outFileNamePrefix 03.Align/$id #^N

samtools index -@ 8 03.Align/${id}\Aligned.sortedByCoord.out.bam  #^N


## counting 
tail name -n 14|while read id;do
#bsub -q normal_1week -n 1 -e $id.err -o $id.log -J $id \
htseq-count -f bam -t exon -i gene_name --minaqual=5 03.Align/$id\Aligned.sortedByCoord.out.bam \
~/genome/mm10/mm10.ucsc.gtf> 04.htseq/$id.count


cut -f 1,7 05.featureCounts/$id.featureCounts.txt|grep -v '^#'>05.featureCounts/$id\_RNA_count.tsv
#merge 
python ~/script/RNAseq/merge_metaphlan_tables.py 05.featureCounts/*NA_count.tsv > 05.featureCounts/merge.all_rna_count.txt


## star pipline
cat name|while read id;do
bsub -q normal_1week -n 8 -e $id.err -o $id.log -J $id \
samtools sort -@ 8 -o 04.Align/$id.sorted.bam 04.Align/$id\Aligned.sortedByCoord.out.bam #^N
samtools index -@ 8 04.Align/$id.sorted.bam #^N



cat name|while read id;do
bsub -q normal -n 20 -e $id.err -o $id.log -J $id \
cufflinks -p 20 -o 05.cufflinks/$id -g /data/gpfs01/wmo/genome/mm10/mm10.ucsc.gtf -u 03.Alignment/$id.sorted.bam

#merge 
python ~/script/RNAseq/merge_metaphlan_tables.py 04.htseq/*.count > 04.htseq/merge.all_ribo.txt



cat name|while read id;do
cp /data/gpfs02/wmo/work/RNA-seq/LHD/LHD_paper/00.rawdata/$id/*sra /data/gpfs02/wmo/work/RNA-seq/LHD/TEST/01.raw
done

scp wmo@hpc.xmu.edu.cn:/data/gpfs02/wmo/work/RNA-seq/ZBX_E18/03.Align/*bam* E18_MIA_bam_v2/

