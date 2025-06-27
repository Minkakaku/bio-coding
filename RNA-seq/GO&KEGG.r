## lw p53
library(AnnotationDbi)
library(org.Mm.eg.db)
library(DOSE)
library(clusterProfiler)
library(magrittr)
library(KEGG.db)
# for R's pipe operator %>%
# clusterProfiler

directory <- "pathtoyourdata"
setwd(directory)
list.files()

d1 <- read.csv("pathtoyourdata", header = TRUE, stringsAsFactors = FALSE)
up <- d1[which(d1$logFC > 1), ]
down <- d1[which(d1$logFC < -1), ]
gene_SYMBOL <- up[, 1]

# Go enrichment_MF
ego <- enrichGO(
  gene = gene_SYMBOL, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
  ont = "ALL", pvalueCutoff = 0.1
)
ego <- ego[order(ego@result$pvalue), ]
write.csv(ego, file = "GO_p35_down.csv")

## pathway
gene_ENTREZID <- unlist(bitr(gene_SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db))
kk <- enrichKEGG(
  gene = gene_ENTREZID,
  organism = "mmu",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, use_internal_data = T
)
kk1 <- kk@result %>% as.data.frame()
SYMBOL <- unlist(lapply(X = kk1$geneID, FUN = function(x) {
  return(paste(unlist(bitr(strsplit(x, split = "/")[[1]],
    fromType = "ENTREZID",
    toType = "SYMBOL",
    OrgDb = org.Mm.eg.db
  )[2]),
  collapse = "/"
  ))
}))
kk1$symbol <- SYMBOL
write.csv(kk1, file = "pathway_p43_down.csv")
