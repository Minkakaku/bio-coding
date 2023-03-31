# repeatM
# Dseq2
rm(list = ls())
library(AnnotationDbi)
library(DESeq2)
library(dplyr)
path <- "~/PJ/YXX/ph40-ph24-aav-lnc/repeatM"
setwd(path)
list.files()

ph24_ko1 <- read.delim("PH24h_KO_1_1_val_1_U_FC_tpm.genloc", header = FALSE)
ph24_ko1 <- unite(ph24_ko1, "ID", c("V1", "V2", "V3", "V8"), remove = TRUE)
ph24_ko1 <- unique(ph24_ko1[, c(1, 4, 10)])
colnames(ph24_ko1) <- c("ID", "ph24_ko1", "geneid")

ph24_ko2 <- read.delim("PH24h_KO_2_1_val_1_U_FC_tpm.genloc", header = FALSE)
ph24_ko2 <- unite(ph24_ko2, "ID", c("V1", "V2", "V3", "V8"), remove = TRUE)
ph24_ko2 <- unique(ph24_ko2[, c(1, 4, 10)])
colnames(ph24_ko2) <- c("ID", "ph24_ko2", "geneid")

ph24_wt1 <- read.delim("PH24h_WT_1_1_val_1_U_FC_tpm.genloc", header = FALSE)
ph24_wt1 <- unite(ph24_wt1, "ID", c("V1", "V2", "V3", "V8"), remove = TRUE)
ph24_wt1 <- unique(ph24_wt1[, c(1, 4, 10)])
colnames(ph24_wt1) <- c("ID", "ph24_wt1", "geneid")

ph24_wt2 <- read.delim("PH24h_WT_2_1_val_1_U_FC_tpm.genloc", header = FALSE)
ph24_wt2 <- unite(ph24_wt2, "ID", c("V1", "V2", "V3", "V8"), remove = TRUE)
ph24_wt2 <- unique(ph24_wt2[, c(1, 4, 10)])
colnames(ph24_wt2) <- c("ID", "ph24_wt2", "geneid")

ph40_ko1 <- read.delim("PH40h_KO_1_1_val_1_U_FC_tpm.genloc", header = FALSE)
ph40_ko1 <- unite(ph40_ko1, "ID", c("V1", "V2", "V3", "V8"), remove = TRUE)
ph40_ko1 <- unique(ph40_ko1[, c(1, 4, 10)])
colnames(ph40_ko1) <- c("ID", "ph40_ko1", "geneid")

ph40_ko2 <- read.delim("PH40h_KO_2_1_val_1_U_FC_tpm.genloc", header = FALSE)
ph40_ko2 <- unite(ph40_ko2, "ID", c("V1", "V2", "V3", "V8"), remove = TRUE)
ph40_ko2 <- unique(ph40_ko2[, c(1, 4, 10)])
colnames(ph40_ko2) <- c("ID", "ph40_ko2", "geneid")

ph40_wt1 <- read.delim("PH40h_WT_1_1_val_1_U_FC_tpm.genloc", header = FALSE)
ph40_wt1 <- unite(ph40_wt1, "ID", c("V1", "V2", "V3", "V8"), remove = TRUE)
ph40_wt1 <- unique(ph40_wt1[, c(1, 4, 10)])
colnames(ph40_wt1) <- c("ID", "ph40_wt1", "geneid")

ph40_wt2 <- read.delim("PH40h_WT_2_1_val_1_U_FC_tpm.genloc", header = FALSE)
ph40_wt2 <- unite(ph40_wt2, "ID", c("V1", "V2", "V3", "V8"), remove = TRUE)
ph40_wt2 <- unique(ph40_wt2[, c(1, 4, 10)])
colnames(ph40_wt2) <- c("ID", "ph40_wt2", "geneid")

f <- list(ph24_ko1, ph24_ko2, ph24_wt1, ph24_wt2, ph40_ko1, ph40_ko2, ph40_wt1, ph40_wt2)
Mergefile <- Reduce(function(x, y) {
  merge(x, y,
    by = "ID", all = TRUE, no.dups = TRUE
  )
}, f)

colnames(Mergefile) <- c(
  "ID", "ph24_ko1", "geneid.x", "ph24_ko2",
  "geneid.1", "ph24_wt1", "geneid.2", "ph24_wt2", "geneid.3", "ph40_ko1",
  "geneid.4", "ph40_ko2", "geneid.5", "ph40_wt1", "geneid.6", "ph40_wt2", "geneid.7"
)
Mergefile2 <- unite(Mergefile, "geneid", c(3, 5, 7, 9, 11, 13, 15, 17), remove = TRUE, sep = ",")
Mergefile2[is.na(Mergefile2)] <- 0
write.csv(Mergefile2, "merge.counts.csv")
Mergefile <- Mergefile2
head(Mergefile)
rownames(Mergefile) <- Mergefile[, 1]
Mergefile <- Mergefile[, -1]
head(Mergefile)
Mergefile <- Mergefile[rowMeans(Mergefile[c(1, 3, 4, 5, 6, 7, 8, 9)]) > 1, ]

countdata <- Mergefile
colnames(countdata)
countdata <- countdata[, c(1, 3, 4, 5, 6, 7, 8, 9)]

condition <- factor(rep(c("ph24ko", "ph24wt", "ph40ko", "ph40wt"), each = 2))

coldata <- data.frame(row.names = colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = coldata,
  design = ~condition
)
dds <- DESeq(dds)
coldata
res <- results(dds, contrast = c("condition", "ph40ko", "ph40wt"))


resordered <- res[order(res$pvalue), ]
# sum(resordered$padj < 0.1, na.rm = TRUE) # 有多少padj小于0.1的
# diff_gene <- subset(resordered, padj < 0.1 & abs(log2FoldChange) > 1)

res <- as.data.frame(resordered)
res$id <- rownames(res)
Mergefile$id <- rownames(Mergefile)
res1 <- merge(res, Mergefile, by = "id")
res1 <- res1[order(res1$pvalue), ]
write.csv(res1, file = "ph40ko|wt.csv", row.names = FALSE)

# gen
genefile = read.delim("/home/hongfan/PJ/YXX/ph40-ph24-aav-lnc/gen/cut/Mergefile")
genefile = genefile[rowMeans(genefile[,c(1:7,9)])>0,]
genefile = unique(genefile)
rownames(genefile) <- genefile[, 8]
genefile <- genefile[, -8]
head(genefile)

countdata <- genefile
colnames(countdata)
condition <- factor(rep(c("ph24ko", "ph24wt", "ph40ko", "ph40wt"), each = 2))

coldata <- data.frame(row.names = colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = coldata,
  design = ~condition
)
dds <- DESeq(dds)
coldata
res <- results(dds, contrast = c("condition", "ph40ko", "ph40wt"))


resordered <- res[order(res$pvalue), ]
# sum(resordered$padj < 0.1, na.rm = TRUE) # 有多少padj小于0.1的
# diff_gene <- subset(resordered, padj < 0.1 & abs(log2FoldChange) > 1)

res <- as.data.frame(resordered)
res$id <- rownames(res)
countdata$id <- rownames(countdata)
res2 <- merge(res, countdata, by = "id")
res2 <- res2[order(res2$pvalue), ]
write.csv(res2, file = "ph40ko|wt|gen.csv", row.names = FALSE)

res1 <- as.data.frame(res1)
res1[,9] <- str_replace(res1[,9], "NA,", "")
res1 = as.data.frame(res1)
res1 = separate(res1,"geneid",c("geneid"),remove = T,extra = 'drop')
res2 = res2[,c(1,3,8,9)]
res3 = merge(res1,res2,by.x="geneid",by.y="id",all.x=T)
write.csv(res3,"repeatMwithgenlfc.csv")

res4 = res3[res3$log2FoldChange.x>1&res3$log2FoldChange.y>1,]
res4id = unique(res4$geneid)
library(clusterProfiler)
library(org.Mm.eg.db)
res4id <- bitr(res4id,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Mm.eg.db)
write.csv(as.data.frame(res4id),"res4id.csv")
ego_all = enrichGO(res4id$ENTREZID,
                   org.Mm.eg.db,
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,readable=T)
library(ggplot2)
write.csv(as.data.frame(ego_all),"ego.csv")
barplot(ego_all,showCategory = 20)
ggsave("ego.pdf",height=20, width = 12)













# plot
# http://www.ehbio.com/ImageGP/index.php/Home/Index/Volcanoplot.html
suppressPackageStartupMessages({
  library(ggpubr)
  library(ggrepel)
  library(gg.gap)
  library(ggplot2)
  library(dplyr)
  library(ggsci)
  library(scales)
  library(tidyr)
})
color3 <- pal_npg()(10)
show_col(color3)
color3
results$id <- rownames(results)
res1 <- results[, c(3, 5, 6)]
colnames(res1) <- c("log2FC", "Pvalue", "gene_id")
res1 <- res1[, c("gene_id", "log2FC", "Pvalue")]
res2 <- res1[complete.cases(res1), ]
head(res2)

res2$change <- as.factor(ifelse(res2$Pvalue < 0.5 & abs(res2$log2FC) > 1,
  ifelse(res2$log2FC > 1, "UP", "DOWN"), "NOT"
))
head(res2)
res2$label <- ifelse(res2$Pvalue < 0.5 & abs(res2$log2FC) > 2, res2$gene_id, "")
head(res2)
# res3$geneid <- res2$gene_id
res2$minuslog10Pvalue <- -log10(res2$Pvalue)
count <- table(res2$change)
a <- ggscatter(res2,
  y = "minuslog10Pvalue",
  x = "log2FC",
  ylab = "-log10(Pvalue)",
  color = "change",
  size = 2,
  # label = "label",
  # font.label = c(12, "#000000"),
  # repel = FALSE,
  palette = c("#4DBBD5FF", "#999999", "#E64B35FF"),
  alpha = 0.6
)
a
a2 <- a + annotate("label",
  x = c(-5, 5), y = 18,
  label = count[c(1, 3)], col = c("#4DBBD5FF", "#E64B35FF")
) +
  annotate("text", x = c(-5, 5), y = 20.5, label = c("DOWN", "UP")) +
  geom_text_repel(aes(label = label),
    size = 4, color = "black",
    point.padding = 0.3,
    fontface = "bold",
    hjust = 0.5,
    arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
    segment.size = 0.1, segment.alpha = 0.8, nudge_y = 0.7,
    max.overlaps = 10000000 # , xlim = c(0, 3.5)
  )
a3 <- a2 + ylab("-log10 (pvalue)") + # 修改y轴名称
  xlab("log2 (FoldChange)") + # 修改x轴名称
  geom_vline(xintercept = c(-1, 1), lty = 3, col = "black", lwd = 0.5) + # 添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.5), lty = 3, col = "black", lwd = 0.5) +
  theme(
    legend.title = element_blank() # 不显示图例标题
  )
a3
ggsave("269wt270mut.pdf", width = 10, height = 10)
