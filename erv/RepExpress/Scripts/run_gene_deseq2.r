rm(list = ls())

library("DESeq2")

control <- "wt"
treatment <- "ko"
countData <- read.csv("gene_count_matrix.csv", header = TRUE, row.names = 1)
countData <- countData[rowMeans(countData) > 1, ]
head(countData)
condition <- factor(c(rep(treatment, 3), rep(control, 3)))
colData <- data.frame(row.names = colnames(countData), condition)
colData
dds <- DESeqDataSetFromMatrix(
    countData = countData,
    colData = colData, design = ~condition
)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dispersions(dds) <- mcols(dds)$dispGeneEst
dds <- nbinomWaldTest(dds)
res <- results(dds, contrast = c("condition", treatment, control))
res_sorted <- res[order(res$padj), ]
write.table(res_sorted,
    file = "gene_deseq2.csv",
    sep = ",", quote = FALSE, col.names = NA, row.names = TRUE, append = TRUE
)
