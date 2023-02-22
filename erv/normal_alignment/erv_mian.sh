# bin/bash
find . -name '*_2.fq.gz' -exec basename {} _2.fq.gz \; > name
# mkdir 00.rawdata 01.qc 02.bam 03.counts ;
# prefetch --option-file name -O 00.rawdata ;

# conda activate py39

for i in `cat name` ; do
export i
bsub -q normal_1day -n 24 -e /data/gpfs02/wmo/work/HF/00.err_log/test.err -o /data/gpfs02/wmo/work/HF/00.err_log/test.log -J 1207test \
bash PC.sh
done
cd 03.counts
Rscript /data/gpfs02/wmo/work/HF/01.scripts/tools/hisat2_merge.r