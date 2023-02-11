##RNASeq with ERV pipeline
## Author: JMY
## Data: 2020-8-9
## Data was from CSX (bam files)

scp -r bam/ wmo@hpc.xmu.edu.cn:/data/gpfs02/wmo/work/repeat-seq/CSX/E18


cd /data/gpfs02/wmo/work/repeat-seq/CSX/E18

# counts
## 参考组有点问题，诺和的用v3，其他的用原版
cat name| while read id 
do
bsub -q normal_1week -n 12 -e $id.err -o $id.log -J $id \
featureCounts -T 12 -t exon -g gene_id -s 1 -a /data/gpfs01/wmo/genome/mm10/rep/ucsc_mm10.rep_v3.gtf \
-o 02.count/$id.txt 01.bam/$id.bam
done


# cut
cat name| while read id 
do
cut -f 1,7 02.count/$id.txt> 02.count/$id.all_TE_count.txt
done

#merge 
python ~/script/RNAseq/merge_metaphlan_tables.py 02.count/*.all_TE_count.txt > 02.count/merge_TE.counts.txt


