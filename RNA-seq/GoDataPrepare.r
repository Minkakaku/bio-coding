rm(list = ls())
# hongfan
suppressPackageStartupMessages({
    library(plyr)
    library(stringr)
    library(ggplot2)
    library(grid)
    library(clusterProfiler)
    library(AnnotationDbi)
    library(org.Mm.eg.db)
})

dir <- "pathtoyourdata"
setwd(dir)
files <- list.files(pattern = "csv")

# circle process
for (i in files) {
    data <- read.csv(i)
    data.name <- unlist(strsplit(i, ".", fixed = TRUE))[1]
    ego <- enrichGO(
        gene = data$gene_id, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",
        ont = "CC", pvalueCutoff = 0.1
    )
    go_data <- ego[order(as.data.frame(ego)$pvalue), ]
    SYMBOL <- unlist(lapply(X = go_data$geneID, FUN = function(x) {
        return(paste(unlist(bitr(strsplit(x, split = "/")[[1]],
            fromType = "ENSEMBL",
            toType = "SYMBOL",
            OrgDb = org.Mm.eg.db
        )[2]),
        collapse = "/"
        ))
    }))
    go_data$symbol <- SYMBOL
    write.csv(go_data, paste0(data.name, "CC.go.csv"))
    go_data$GeneRatio <- mixedToFloat(go_data$GeneRatio)

    go_data$set <- as.numeric(set)
    genenumber <- as.numeric(go_data$Count)
    # df$GO.Term <- paste0(str_to_upper(str_sub(df$GO.Term,1,1)), str_sub(df$GO.Term,2,str_length(df$GO.Term)))
}
