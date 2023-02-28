library("DESeq2")
countData <- read.csv("DESeq2_count_matrix.txt", header=TRUE, sep="\t")
metaData <- read.csv("DESeq2_metadata.txt", header=TRUE, sep="\t")
dds <- DESeqDataSetFromMatrix(countData=countData,colData=metaData,design=~descr,tidy=TRUE)
dds <- DESeq(dds)
res <- results(dds, contrast = c("descr", "treatment", "control"))
res_sorted <- res[order(res$padj),]
write.table(res_sorted,file="../05.deseq2_out/EWAT_deseq2_express.txt",sep=",",quote=FALSE,col.names=NA,row.names=TRUE,append=TRUE)
q()
