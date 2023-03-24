rm(list = ls())
library(AnnotationDbi)
library(DESeq2)
library(dplyr)


path <- "/home/hongfan/Punchure/projects/YHF_org/03.counts"
setwd(path)
list.files()

countdata <- read.csv("merge.counts.csv", header = TRUE, row.names = 2)
countdata <- countdata[, -1]
countdata <- countdata[rowMeans(countdata) > 1, ]
c("24", "48", "72", "Mko_WEo", "Mko_WEe")

countdata <- countdata[, grep("48", colnames(countdata))]

treat <- as.data.frame(strsplit(colnames(countdata), "h"))[1, ]
treat <- unique(as.character(treat))
condition <- factor(rep(treat, each = 2))
# time <- factor(c(
#     rep(rep(c("24", "72"), each = 2), 4),
#     rep(rep(c("24", "48", "72"), each = 2), 2),
#     rep(rep(c("24", "72"), each = 2), 2)
# ))
# coldata <- data.frame(row.names = colnames(countdata), condition, time)
coldata <- data.frame(row.names = colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(
    countData = countdata,
    colData = coldata,
    design = ~condition
)
dds <- DESeq(dds)
coldata
# KO 在前
res_24 <- results(dds, contrast = c("condition", "Mko_WEo_24", "Mko_WEe_24"))
res_72 <- results(dds, contrast = c("condition", "Mko_WEo_72", "Mko_WEe_72"))
res_48 <- results(dds, contrast = c("condition", "Mko_WEo_48", "Mko_WEe_48"))

resordered <- res_48[order(res_48$pvalue), ]
# sum(resordered$padj < 0.1, na.rm = TRUE) # 有多少padj小于0.1的
# diff_gene <- subset(resordered, padj < 0.1 & abs(log2FoldChange) > 1)

res <- as.data.frame(resordered)
res$id <- rownames(res)
countdata <- countdata[, grep("Mko_WE", colnames(countdata))]
countdata <- countdata[, grep("48", colnames(countdata))]
countdata$id <- rownames(countdata)
res1 <- merge(res, countdata, by = "id")
res1 <- res1[order(res1$pvalue), ]
write.csv(res1, file = "res_48.csv", row.names = FALSE)
