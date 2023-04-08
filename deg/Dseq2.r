rm(list = ls())
library(AnnotationDbi)
library(DESeq2)
library(dplyr)
library(tidyr)

path <- "/home/hongfan/PJ/YXX/NASH/repeat"
setwd(path)
list.files()

countdata <- read.delim("merge.txt", header = TRUE)
head(countdata)
countdata <- countdata[, -6]
countdata = unite(countdata, "id",
    c("Chromosome", "start", "end", "RE_uniq_ID"),
    sep = ","
)
head(countdata)

exp = countdata[,c(1,11:22)]
head(exp)
rownames(exp) = exp[, 1]
exp = exp[,-1]
exp <- exp[rowMeans(exp) > 1, ]

coldata = read.csv("SraRunTable.csv")
coldata =coldata[,c("Run",'Phenotype')]
coldata

dds <- DESeqDataSetFromMatrix(
    countData = exp,
    colData = coldata,
    design = ~Phenotype
)
dds <- DESeq(dds)
coldata
# KO 在前
res <- results(dds, contrast = c("Phenotype", "NASH", "chow"))
resordered <- res[order(res$pvalue), ]
# sum(resordered$padj < 0.1, na.rm = TRUE) # 有多少padj小于0.1的
# diff_gene <- subset(resordered, padj < 0.1 & abs(log2FoldChange) > 1)

res <- as.data.frame(resordered)
res$id <- rownames(res)
res1 <- merge(res, countdata, by = "id")
res1 <- res1[order(res1$pvalue), ]
write.csv(res1, file = "NASH,CHOW,repeatM.csv", row.names = FALSE)
