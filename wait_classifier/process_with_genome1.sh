#!/bin/bash
set -e
set -u
set -o pipefail
## create gff file
cd /data/gpfs01/wmo/genome/mm10/
#bsub -q normal -n 5 -e test.err -o test.log -J test RepeatMasker -pa 5 -no_is -species "mus musculus" -html -gff -dir rep/ test.fa
bsub -q normal_1day -n 20 -e Re.err -o Re.log -J Re \
RepeatMasker -pa 20 -no_is -species "mus musculus" -html -gff -dir rep/ mm10.fa
 

sed '1,3d' mm10.fa.out > mm10_d.fa.out

#awk '{print $10,$11}' mm10_d.fa.out > mm10_d_class.txt
#注:需要反复应用awk、sed等shell命令自建符合自己要求的ref database
awk '{if($1>=2000); print $5"\t"$11"\t""exon""\t"$6"\t"$7"\t"$2"\t"$9"\t"".""\t""gene_id ""\""$10"\""}' mm10_d.fa.out > ucsc_mm10.rep.gtf
awk '{if($1>=2000); print $10"\t"$11}' mm10.fa.out > ucsc_mm10.temp  
sort ucsc_mm10.temp > ucsc_mm10.ref.temp 
uniq ucsc_mm10.ref.temp > ucsc_mm10.ref