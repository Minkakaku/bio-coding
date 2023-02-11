require(clusterProfiler)
require(org.Mm.eg.db)
# papare data
d <- "pathtoyourdata"
setwd(d)
df <- read.csv("file", header = TRUE)
keytypes(org.Mm.eg.db)
# methods
eg <- bitr(df$ensembl_id,
    fromType = "ENSEMBL",
    toType = "SYMBOL",
    OrgDb = "org.Mm.eg.db"
)

results <- merge(df, eg,
    by.x = "ensembl_id",
    by.y = "ENSEMBL",
    all.x = T
)

id <- na.omit(results$SYMBOL)
df.entrezID <- mapIds(org.Mm.eg.db,
    id,
    "SYMBOL",
    column = "ENTREZID"
)
