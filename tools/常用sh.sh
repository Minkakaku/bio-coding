#!/usr/bin/bash
grep LINE rmsk.txt | awk -F"\t" '{print $6,$7,$8,$10}' | sed 's/ /\t/g' >LINE.bed
awk 'NR>1 {print $6,$7,$8,$10,"0",$11}' rmsk.txt >mm39.TE.bed

bsub -q normal -n 12 -K -e /data/gpfs02/wmo/work/HF/err_log/1.err -o /data/gpfs02/wmo/work/HF/err_log/1.log -J scTE_build \
    scTE_build -te /data/gpfs02/wmo/genome/mm39/mm39.TE.bed -gene /data/gpfs02/wmo/genome/mm39/genes/refGene.gtf -o /data/gpfs02/wmo/genome/mm39/TE_index
scTE_build -g mm10

bsub -m 'c045 c046' -n 40 -e /data/gpfs02/wmo/work/HF/err_log/2.err -o /data/gpfs02/wmo/work/HF/err_log/2.log -J model_build \
    python /data/gpfs02/wmo/work/HF/scripts/scRNA-seq/model_build.py

bsub -q normal -n 24 -e /data/gpfs02/wmo/work/HF/err_log/1.err -o /data/gpfs02/wmo/work/HF/err_log/1.log -J scTE_build

bsub -q fat4 -n 20 -e 2.err -o 2.log -J model_build \
    python /share/home/zju_qiumj/workbench/tools/scRNA-seq/model_build.py
############################
scTE -i /data/gpfs02/wmo/work/HF/scRNA-seq/HF/20220727scTE/possorted_genome_bam.bam -o /data/gpfs02/wmo/work/HF/scRNA-seq/HF/20220727scTE -x /data/gpfs02/wmo/genome/hg38/TE/hg38.exclusive.idx -CB CB -UMI UB

bamtools split -in /data/gpfs02/wmo/work/HF/scRNA-seq/HF/20220727scTE/possorted_genome_bam.bam -tag CB >filted.bam

bsub -q normal_1week -n 24 -e /data/gpfs02/wmo/work/HF/err_log/2.err -o /data/gpfs02/wmo/work/HF/err_log/2.log -J R \
    Rscript /data/gpfs02/wmo/work/HF/scripts/scRNA-seq/cpdbprepare.r

bsub -q normal_1week -n 24 -e /data/gpfs02/wmo/work/HF/err_log/aged_S.err -o /data/gpfs02/wmo/work/HF/err_log/aged_S.log -J aged_S \
    cellphonedb method statistical_analysis aged_meta.txt aged_S --counts-data hgnc_symbol --output-path=aged_S_out

bsub -m 'c045 c046 c047 c048' -n 5 -e /data/gpfs02/wmo/work/HF/err_log/p1.err -o /data/gpfs02/wmo/work/HF/err_log/p1.log -J p1 \
    Rscript /data/gpfs02/wmo/work/HF/scripts/scRNA-seq/schep.r

bsub -q normal_1week -K -n 24 -e /data/gpfs02/wmo/work/HF/err_log/p1.err -o /data/gpfs02/wmo/work/HF/err_log/p1.log -J p1 \
    Rscript /data/gpfs02/wmo/work/HF/01.scripts/scRNA-seq/sct_integrate.r
