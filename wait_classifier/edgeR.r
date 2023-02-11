## lhd cm d123 rapa
library("AnnotationDbi")
library("org.Mm.eg.db")
library("edgeR")

directory <- "/home/molab/Share/hf/Project/YHF/04.featureCounts"
setwd(directory)
list.files()

data <- read.table("merge.all_RNA_counts.txt",
    header = T,
    fill = TRUE
)


data1 <- as.data.frame(sapply(data[, c(3, 5)], as.numeric))
data1[is.na(data1)] <- 0
rownames(data1) <- data$ID
group <- c(2, 1)


y <- DGEList(count = data1, group = group) ## mut放在后面

keep <- rowSums(cpm(y) > 1) >= 2 # 至少在两个样本里cpm大于1
y <- y[keep, , keep.lib.sizes = FALSE]

y <- calcNormFactors(y)
y$samples
# plotMDS(y)

# bcv:如果是人类数据, 且实验做的很好,设置为0.4, 如果是遗传上相似的模式物种, 设置为0.1, 如果是技术重复, 那么设置为0.01
y_bcv <- y
bcv <- 0.1
et2 <- exactTest(y_bcv, dispersion = bcv^2)
gene2 <- decideTestsDGE(et2, p.value = 0.05, lfc = 0)
summary(gene2)

res <- et2$table[order(et2$table$PValue), ]
res2 <- data.frame(ID = rownames(res), res[, 1:3])
res3 <- merge(data[, c(1, 2, 6)], res2, by.x = "ID", by.y = "ID")
res <- res3[order(res3$PValue), ]

write.csv(res, file = "5vs3.csv")


anno <- read.csv(file = "anno_TE.csv", header = T, stringsAsFactors = F)

data <- read.csv(file = "diff_rapa_D3mut_D3wt.csv", header = T, stringsAsFactors = F)

all <- merge(data[, 2:7], anno, by.x = "ID", by.y = "V1")
all1 <- all[order(all$PValue), ]

write.csv(all1, file = "diff_rapa_D3mut_D3wt_anno.csv", row.names = F)

list.files()
