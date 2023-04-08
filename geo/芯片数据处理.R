rm(list = ls())
dir = "/home/hongfan/Punchure/projects/CQ_NASHvsHC"
setwd(dir)
list.files()
library(GEOquery)
gse_number <- "GSE89632"
eset <- getGEO(gse_number,
    destdir = "00.rawdata",
    getGPL = TRUE
)
class(eset)
length(eset)
eset <- eset[[1]]
exp <- exprs(eset)
exp[1:4, 1:4]
if (FALSE) {
    exp <- log2(exp + 1)
}
boxplot(exp)
if (FALSE) {
    exp <- limma::normalizeBetweenArrays(exp)
}
boxplot(exp)
pd <- pData(eset)
gpl_number <- eset@annotation
library(tidyverse)
pd.new <- separate(pd,
    col = "title",
    into = c("org", "condition", "number"), sep = "_"
)
group <- pd.new[grep("NASH|HC", pd.new$condition), ]
p <- identical(rownames(group), colnames(exp))
p
if (!p) exp <- exp[, match(rownames(group), colnames(exp))]
name <- rownames(group)
group = group[, "condition"]
group <- factor(group,
    levels = c("HC", "NASH")
)

# library(tinyarray)
# find_anno(gpl_number)
# gpl_number
# ids <- AnnoProbe::idmap("GPL14951")
###################################################################################
if (TRUE) {
    library(GEOquery)
    # 注：soft文件列名不统一，活学活用，有的表格里没有symbol列，也有的GPL平台没有提供注释
    a <- getGEO(gpl_number, destdir = ".")
    b <- a@dataTable@table
    colnames(b)
    ids2 <- b[, c("ID", "Symbol")]
    colnames(ids2) <- c("probe_id", "symbol")
    ids2 <- ids2[ids2$symbol != "" & !str_detect(ids2$symbol, "///"), ]
}


library(limma)
design <- model.matrix(~group)
fit <- lmFit(exp, design)
fit <- eBayes(fit)
deg <- topTable(fit, coef = 2, number = Inf)

pvalue_t <- 0.05
logFC_t <- 1
k1 <- (deg$P.Value < pvalue_t) & (deg$logFC < -logFC_t)
table(k1)
k2 <- (deg$P.Value < pvalue_t) & (deg$logFC > logFC_t)
table(k2)
deg$change <- ifelse(k1, "DOWN", ifelse(k2, "UP", "NOT"))
table(deg$change)
deg <- merge(deg, exp, by = "row.names")
library(dplyr)
deg <- mutate(deg, probe_id = rownames(deg))
ids2 <- ids2[!duplicated(ids2$symbol), ]
deg2 <- inner_join(deg, ids2, by = "probe_id")
# 将deg与ids合并  自动把deg中重复的probe_id给删除
nrow(deg)
