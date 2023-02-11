
##refexpress使用步骤
##来自于文章Stockwell, P. A., Lynch-Sutherland, C. F., Chatterjee, A., Macaulay,
##E. C., & Eccles, M. R. (2021). RepExpress: A novel pipeline for the
##quantification and characterization of transposable element
##expression from RNA-seq data. Current Protocols, 1, e206.
##doi: 10.1002/cpz1.206
#1.准备工作 从UCSC上下载对应生物的全基因组，以及基因组注释文件（GTF），同时从https://genome.ucsc.edu/cgi-bin/hgTables上下载带有REPEATMASKER ANNOTATION的GTF文件

#RM ANNO GTF设置要求如下
#Downloads of TEs from repeatmasker tracks from UCSC are available from: https://genome.ucsc.edu/cgi-bin/hgTables, but some steps are necessary to retrieve the data in an appropriate and complete form.
#assembly: must be set to the required build (hg19)
#group: must be set to 'Repeats'
#output format: must be set to 'all fields from selected table'
#in order to retrieve details that are omitted from the 'GTF - gene trans-fer format (limited)' format.
#output file: should be set to something to avoid having the data sent directly for display by the browser (e.g. hg19_ucsc_repeats.txt)

#或者，也可以（wget http://hgdownload.cse.ucsc.edu/goldenPath/mm39/database/rmsk.txt.gz）
#2. 建立索引（INDEX），每一种生物只需要制作执行一次建立索引工作，将basic_params.sh中的参数设置完全后，执行如下命令。
bsub -q normal_1day_new -n 4 -e TEST1.err -o TEST.log -J TEST \
./generate_combined_gtfs.sh ./basic_params.sh

#3.回帖和计数（注意输入的rawdata必须是fastq格式，如果是fq.gz会报错）(可以用QC后的FQ文件进行)
bsub -q normal_1week -n 4 -e MAP.err -o MAP.log -J MAP \
./map_count_reads.sh ./basic_params.sh ./map_params3.sh

#4.将不同样品的表达矩阵汇总
./compare_express.sh ./basic_params.sh ./combine_params.sh

#5.进行组织富集和差异富集（RETE）(需要准备包含样本文件tpm.genloc路径的genloc.list.txt)
./contrast_express.sh ./basic_params.sh ./combine_params.sh

#6.进行DESEQ2差异分析(需要准备包含样本文件U_FC.tpm路径的tpm.list.txt)
./run_deseq2.sh ./basic_params.sh ./DESeq2_params.sh


 

