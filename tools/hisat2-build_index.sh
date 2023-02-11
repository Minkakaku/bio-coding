bsub -q normal_1day -n 24 -e /data/gpfs02/wmo/work/HF/00.err_log/$i.err -o /data/gpfs02/wmo/work/HF/00.err_log/$i.log -J 1128 \
hisat2-build -f /data/gpfs02/wmo/genome/mm39/mm39.fa -p 24 \
--ss /data/gpfs02/wmo/genome/mm39/mm39.ss \
--exon /data/gpfs02/wmo/genome/mm39/mm39.exon \
/data/gpfs02/wmo/genome/hisat2_index/mm39/mm39