##RNASeq with ERV pipeline
## Author: JMY
## Data: 2021.9.22
## Data was from paper

cd /data/gpfs02/wmo/work/RNA-seq/LHD/LHD_paper

# counts
cat name| while read id ;do
bsub -q normal_1week -n 12 -e $id.err -o $id.log -J $id \
featureCounts -T 12 -t exon -g gene_id -s 1 -a /data/gpfs01/wmo/genome/mm10/rep/ucsc_mm10.rep_v3.gtf \
-o 05.erv_count/$id.txt 03.Alignment/$id.sorted.bam
done

cat name| while read id 
do
cut -f 1,7 $id.txt> $id.all_count.txt
done

#merge 
python ~/script/RNAseq/merge_metaphlan_tables.py 05.erv_count/*count.txt > 05.erv_count/merge.all_TE_counts.txt





