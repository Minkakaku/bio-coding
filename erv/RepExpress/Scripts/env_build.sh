mkdir 00.rawdata 01.qc 02.mapped_out 03.combined_out 04.contrast_out 05.deseq2_out
mv *fq.gz 00.rawdata/
cd 00.rawdata/
find . -name '*_2.fq.gz' -exec basename {} _2.fq.gz \; > name