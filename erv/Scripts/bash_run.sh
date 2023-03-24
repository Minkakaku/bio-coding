cd Scripts

bsub -q normal_1day_new -n 24 -e /data/gpfs02/wmo/HF_work/00.err_log/generate_combined_gtfs.err \
-o /data/gpfs02/wmo/HF_work/00.err_log/generate_combined_gtfs.log -J generate_combined_gtfs \
./generate_combined_gtfs.sh basic_params.sh

for i in `cat name`;do
bsub -q normal_1day_new -n 8 -e /data/gpfs02/wmo/HF_work/00.err_log/$i.err \
-o /data/gpfs02/wmo/HF_work/00.err_log/$i.log -J $i \
trim_galore -q 20 -j 8 --phred33 --stringency 3 --length 20 -e 0.1 --fastqc\
            --paired ../00.fa/$i\_1.fq.gz ../00.fa/$i\_2.fq.gz  \
            --gzip -o ../01.qc
done

for i in `cat name`;do
export i
bsub -q normal_1day_new -n 6 -e /data/gpfs02/wmo/HF_work/$i.err \
-o /data/gpfs02/wmo/HF_work/$i.log -J $i \
./map_count_reads.sh basic_params.sh map_params.sh
done

bsub -q normal_1week -n 24 -e /data/gpfs02/wmo/HF_work/00.err_log/run_deseq2.err \
-o /data/gpfs02/wmo/HF_work/00.err_log/run_deseq2.log -J run_deseq2 \
./run_deseq2.sh basic_params.sh deseq2_params.sh

bsub -q normal_1day_new -n 24 -e /data/gpfs02/wmo/HF_work/00.err_log/samtools.err \
-o /data/gpfs02/wmo/HF_work/00.err_log/samtools.log -J samtools \
samtools index -@ 24 *.bam