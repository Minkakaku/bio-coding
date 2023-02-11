library("AnnotationDbi")
library("org.Mm.eg.db")
library("DESeq2")
library(dplyr)


directory <- "/data/gpfs02/wmo/work/HF/ERV/WM/02.count"
setwd(directory)
list.files()

countData <- read.table("merge_TE.counts.txt",
    header = T,
    stringsAsFactors = F,
    sep = "\t"
)
colData <- read.table("group",
    header = T,
    stringsAsFactors = F,
    sep = "\t"
)

index <- duplicated(countData$ID)
countData <- countData[!index, ]

countData1 <- as.data.frame(sapply(countData[, 2:11], as.numeric))
countData1[is.na(countData1)] <- 0
rownames(countData1) <- countData$ID

a <- c(5:10)
# ctxMH
dds <- DESeqDataSetFromMatrix(
    countData = countData1[, a],
    colData = colData[a, ],
    design = ~Group
)
head(assay(dds))
ddsTC <- DESeq(dds)
sizeFactors(ddsTC)

# 要把wt放后面，因为是前面除以后面来比对的
res <- results(ddsTC, contrast = c("Group", "HF", "SHAM_HF"))
res <- as.data.frame(res[order(res$pvalue), ])
summary(res)
table(res$pvalue < 0.05)
res$ID <- rownames(res)
countData1$ID <- rownames(countData1)
res1 <- merge(res, countData1[, c(a, 11)], by.x = "ID", by.y = "ID")
res1 <- res1[order(res1$pvalue), ]

write.csv(res1, file = "diff_RNA_HF_SHAMHF.csv", row.names = F)

## anno
directory <- "D:/Users/24432/Desktop/生信/RNA&ERV/LHD/paper/TE"
setwd(directory)
list.files()

anno <- read.csv("anno_TE.csv", header = T, stringsAsFactors = F)
data <- read.csv("diff_TE_TAC_ctrl.csv", header = T, stringsAsFactors = F)

d <- merge(data, anno, by.x = "ID", by.y = "V1")

write.csv(d, file = "diff_TE_TAC_ctrl_withanno.csv", row.names = F)