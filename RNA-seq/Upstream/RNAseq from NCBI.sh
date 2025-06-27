## RNA-seq for paper
## data from http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155972

## 1. download srr file
## create new folders
mkdir 00.rawdata/ 01.fasta/ 02.QC/ 03.Alignment/ 04.featureCounts/
prefetch --option-file name -O 00.rawdata/

## translate srr file to fastq file

cat SRR_Acc_List |while read id; do
bsub -q normal -n 1 -e $id.err -o $id.log -J $id   \
 fastq-dump 00.rawdata/$id/$id.sra --gzip --split-e --defline-seq '@$ac-$si/$ri' --defline-qual '+' \
 -O 01.fastq/$id
done

## QC
# 单
cat SRR_Acc_List |while read id;do
bsub -q normal_1day_new -n 12 -e $id.err -o $id.log -J $id \
fastp -w 20 -i 01.fastq/$id/$id.fastq.gz -o 02.QC/${id}.clean.fq.gz \
-j 02.QC/$id.fastp.json -h 02.QC/$id.fastp.html
done

# 双
cat SRR_Acc_List |while read id;do
bsub -q normal_1day_new -n 12 -e $id.err -o $id.log -J $id \
fastp -w 20 -i 01.fastq/$id/$id\_1.fastq.gz -I 01.fastq/$id/$id\_2.fastq.gz \
-o 02.QC/${id}\_1.clean.fq.gz -O 02.QC/${id}\_2.clean.fq.gz \
-j 02.QC/$id.fastp.json -h 02.QC/$id.fastp.html
done

## Align to Genome, sort and index bam file
# 名字短的75bp，名字长的100bp（索引已经建立，此步可不做）

bsub -q normal_1day_new -n 24 -e ind.err -o ind.log -J ind \
STAR --runMode genomeGenerate --genomeDir /home/molab/Share/hf/seqref/STAR_index_GRCh38.p13_300bp \
--genomeFastaFiles /home/molab/Share/hf/GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna --sjdbGTFfile /home/molab/Share/hf/GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf \
--runThreadN 10 --sjdbOverhang 299

bsub -q normal_1day_new -n 10 -e ind.err -o ind.log -J ind \
STAR --runMode genomeGenerate --genomeDir /home/molab/Share/hf/seqref/STAR_index_mm10_55bp \
--genomeFastaFiles /home/molab/Share/hf/mm10.fa --sjdbGTFfile /home/molab/Share/hf/mm10.ucsc.gtf \
--runThreadN 10 --sjdbOverhang 49

# 比对
cat SRR_Acc_List|while read id;do
bsub -q normal_1day_new -n 24 -e $id.err -o $id.log -J $id \
STAR --runThreadN 16 --readFilesCommand zcat --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate \
--genomeDir /home/molab/Share/hf/seqref/STAR_index_mm10_45bp --readFilesIn 02.QC/${id}.clean.fq.gz \
--outFileNamePrefix 03.Alignment/$id
done

# 排序
cat SRR_Acc_List|while read id;do
bsub -q normal_1week -n 8 -e $id.err -o $id.log -J $id \
samtools sort -@ 20 -o 03.Alignment/$id.sorted.bam 03.Alignment/$id\Aligned.sortedByCoord.out.bam
done

# 建立索引，当用XX查看时需要，不是必需
cat SRR_Acc_List|while read id;do
samtools index 03.Alignment/$id.sorted.bam
done


## counting 
cat SRR_Acc_List|while read id;do
bsub -q normal -n 20 -e $id.err -o $id.log -J $id \
featureCounts -T 30 -p -t exon -g gene_id -a /home/molab/Share/hf/mm10.ucsc.gtf -o 04.featureCounts/$id.featureCounts.txt 03.Alignment/$id.sorted.bam
done

# 取ID和counts
cat SRR_Acc_List|while read id;do
cut -f 1,7 04.Counts/$id.Counts.txt|grep -v '^#'>04.Counts/$id\.Counts.filter.tsv
done

#merge 
python /home/molab/Share/hf/Scripts/merge_metaphlan_tables.py 04.Counts/*.tsv > 04.Counts/merge.all_RNA_counts.txt
join 04.featureCounts/GSE73723merge.all_RNA_counts.txt 04.featureCounts/SRR2557086.featureCounts.filter.tsv > 04.featureCounts/GSE73723merge.all_RNA_counts.txt