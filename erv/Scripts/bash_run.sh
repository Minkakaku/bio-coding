cd Scripts

bsub -q normal_1day_new -n 24 -e generate_combined_gtfs.err \
-o generate_combined_gtfs.log -J generate_combined_gtfs \
./generate_combined_gtfs.sh basic_params.sh

for i in `cat name`;do
bsub -q normal_1day_new -n 8 -e $i.err \
-o $i.log -J $i \
trim_galore -q 20 -j 8 --phred33 --stringency 3 --length 20 -e 0.1 --fastqc\
            --paired ../00.fa/$i\_1.fq.gz ../00.fa/$i\_2.fq.gz  \
            --gzip -o ../01.qc
done

for i in `cat name`;do
export i
bsub -q normal_1day_new -n 24 -e $i.err -o $i.log -J $i \
./map_count_reads.sh basic_params.sh map_params.sh
done

bsub -q normal_1week -n 24 -e run_deseq2.err \
-o run_deseq2.log -J run_deseq2 \
./run_deseq2.sh basic_params.sh deseq2_params.sh

bsub -q normal_1day_new -n 24 -e samtools.err \
-o samtools.log -J samtools \
samtools index -@ 24 *.bam