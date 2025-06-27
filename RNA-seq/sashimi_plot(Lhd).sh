##RNASeq pipeline
## Author: JMY
## Data: 2020-5-30
## Data was from Lhd 


cd /data/gpfs01/wmo/work/RNA-seq/Lhd_crypt
## MISO for sashimi plot
# installed by conda: conda install -c bioconda misopy

#gff3 was downloaded from https://miso.readthedocs.io/en/fastmiso/annotation.html
# build index
mkdir -p 02.sashimi_plot/01.miso 02.sashimi_plot/02.summaries 02.sashimi_plot/03.comparison_output \
02.sashimi_plot/04.filtered 02.sashimi_plot/05.sashimi
## create index for gff3
# gff3 classes :1.Skipped exons (SE)  2.Alternative 3’/5’ splice sites (A3SS, A5SS)  3.Mutually exclusive exons (MXE)  4.Tandem 3’ UTRs (TandemUTR)  5.Retained introns (RI)  6.Alternative first exons (AFE)  7.Alternative last exons (ALE) 
#bsub -q normal_1day -n 1 -e index.err -o index.log -J index index_gff --index miso_annotations_mm10_v2/SE.mm10.gff3 ~/Index/MISO_SE_mm10
#index_gff --index miso_annotations_mm10_v2/A5SS.mm10.gff3 ~/Index/MISO_A5SS_mm10
index_gff --index miso_annotations_mm10_v2/RI.mm10.gff3 ~/Index/MISO_RI_mm10
# index_gff --index miso_annotations_mm10_v2/A3SS.mm10.gff3 ~/Index/MISO_A3SS_mm10


## run miso
cat md.txt |while read id
do
bsub -q normal_1day -n 1 -e ${id%%.*}.err -o ${id%%.*}.log -J ${id%%.*} miso --run ~/Index/MISO_RI_mm10/ 01.bam_data/$id \
--output-dir 02.sashimi_plot/01.miso/${id%%.*} --read-len 50 --settings-filename 02.sashimi_plot/miso_settings.txt
## summarize
# bsub -q normal_1day -n 1 -e $id.err -o $id.log -J $id summarize_miso --summarize-samples 05.sashimi_plot/01.miso/$id/ 05.sashimi_plot/02.summaries/$id/
## comparation
# compare_miso --compare-samples control/ knockdown/ comparisons/
done

## Make pairwise comparisons between samples to detect differentially expressed isoforms/events 
# bsub -q normal_1day -n 1 -e summ.err -o summ.log -J summ summarize_miso --summarize-samples 02.sashimi_plot/01.miso/${id%%.*} 02.sashimi_plot/02.summaries/

##确定位点，index里的gff文件
cd ~/Index/MISO_RI_mm10
grep -n "ENSMUSG00000027514" genes.gff

## Number of mapped reads in each sample
cat md.txt |while read id
do
bsub -q normal_1day -n 1 -e ${id%%.*}.err -o ${id%%.*}.log -J ${id%%.*} samtools view -c -F 260 01.bam_data/$id
done

# draw sashimi plot
cd ~/work/RNA-seq/Lhd_crypt

bsub -q normal -n 1 -e plot.err -o plot.log -J plot sashimi_plot --plot-event "chr2:173212266-173212099:-@chr2:173210733-173210536:-" ~/Index/MISO_RI_mm10/ 02.sashimi_plot/sashimi_plot_settings.txt --output-dir 02.sashimi_plot/05.sashimi/



