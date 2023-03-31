mkdir 00.fa 01.qc 02.map
mv *fq.gz 00.fa/
cd 00.fa/
find . -name '*_2.fq.gz' -exec basename {} _2.fq.gz \; > name