library(dplyr)
d <- "~/Share/hf/Project/ZHR/counts"
setwd(d)
list.files()

counts <- read.table("GSE52564merge.all_RNA_counts.txt",
    stringsAsFactors = FALSE,
    header = TRUE,
    fill = TRUE,
    row.names = 1
)
len <- read.delim("GSE52564merge.all_RNA_lengths.txt",
    stringsAsFactors = FALSE,
    header = TRUE,
    fill = TRUE
)
lens_name <- len$ID
lengths <- as.data.frame(lapply( # nolint
    len, # nolint
    as.numeric
))
lens_num <- select(lengths, 2)
lens_num <- unlist(lens_num)

# ng <- intersect(rownames(counts), len$ID)
# length(ng)
# lens <- len[match(ng, len$ID), 2]
names(lens_num) <- lens_name
head(lens_num)

counts[, 1:ncol(counts)] <- as.data.frame(lapply( # nolint
    counts[, 1:ncol(counts)], # nolint
    as.numeric
))

total_counts <- colSums(counts[, 1:ncol(counts)], na.rm = TRUE, dims = 1) # nolint
total_counts[1]

nm_fpkm <- t(do.call(
    rbind,
    lapply(
        1:length(total_counts), # nolint
        function(i) {
            10^9 * counts[, i] / lens_num / total_counts[i] # nolint
        }
    )
))
counts_colname <- colnames(counts)
colnames(nm_fpkm) <- counts_colname
write.csv(nm_fpkm, file = "gse52564.csv")

# 文档结束————————————————————————————————————————————————————————————————————————————————————————————————————

counts <- read.delim("GSE73721hmmerge.all_RNA_counts.txt",
    header = TRUE,
    fill = TRUE,
    row.names = 1
)
len <- read.delim("GSE73721hmmerge.all_RNA_lengths.txt",
    header = TRUE,
    fill = TRUE
)
lens_name <- len$ID
lengths <- as.data.frame(lapply( # nolint
    len, # nolint
    as.numeric
))
lens_num <- select(lengths, 2)
lens_num <- unlist(lens_num)

# ng <- intersect(rownames(counts), len$ID)
# length(ng)
# lens <- len[match(ng, len$ID), 2]
names(lens_num) <- lens_name
head(lens_num)

counts[, 1:ncol(counts)] <- as.data.frame(lapply( # nolint
    counts[, 1:ncol(counts)], # nolint
    as.numeric
))

total_counts <- colSums(counts[, 1:ncol(counts)], na.rm = TRUE, dims = 1) # nolint

total_counts[1]

nm_fpkm <- t(do.call(
    rbind,
    lapply(
        1:length(total_counts), # nolint
        function(i) {
            10^9 * counts[, i] / lens_num / total_counts[i] # nolint
        }
    )
))
counts_colname <- colnames(counts)
colnames(nm_fpkm) <- counts_colname
write.csv(nm_fpkm, file = "gse73721hm.csv")

# 文档结束————————————————————————————————————————————————————————————————————————————————————————————————————



counts <- read.delim("GSE73721mmmerge.all_RNA_counts.txt",
    header = TRUE,
    fill = TRUE,
    row.names = 1
)
len <- read.delim("GSE73721mmmerge.all_RNA_lengths.txt",
    header = TRUE,
    fill = TRUE
)
lens_name <- len$ID
lengths <- as.data.frame(lapply( # nolint
    len, # nolint
    as.numeric
))
lens_num <- select(lengths, 2)
lens_num <- unlist(lens_num)

# ng <- intersect(rownames(counts), len$ID)
# length(ng)
# lens <- len[match(ng, len$ID), 2]
names(lens_num) <- lens_name
head(lens_num)

counts[, 1:ncol(counts)] <- as.data.frame(lapply( # nolint
    counts[, 1:ncol(counts)], # nolint
    as.numeric
))

total_counts <- colSums(counts[, 1:ncol(counts)], na.rm = TRUE, dims = 1) # nolint

total_counts[1]

nm_fpkm <- t(do.call(
    rbind,
    lapply(
        1:length(total_counts), # nolint
        function(i) {
            10^9 * counts[, i] / lens_num / total_counts[i] # nolint
        }
    )
))
counts_colname <- colnames(counts)
colnames(nm_fpkm) <- counts_colname
write.csv(nm_fpkm, file = "GSE73721mm.csv")

# 文档结束————————————————————————————————————————————————————————————————————————————————————————————————————


counts <- read.delim("GSE75246merge.all_RNA_counts.txt",
    header = TRUE,
    fill = TRUE,
    row.names = 1
)
len <- read.delim("GSE75246merge.all_RNA_lengths.txt",
    header = TRUE,
    fill = TRUE
)
lens_name <- len$ID
lengths <- as.data.frame(lapply( # nolint
    len, # nolint
    as.numeric
))
lens_num <- select(lengths, 2)
lens_num <- unlist(lens_num)

# ng <- intersect(rownames(counts), len$ID)
# length(ng)
# lens <- len[match(ng, len$ID), 2]
names(lens_num) <- lens_name
head(lens_num)

counts[, 1:ncol(counts)] <- as.data.frame(lapply( # nolint
    counts[, 1:ncol(counts)], # nolint
    as.numeric
))

total_counts <- colSums(counts[, 1:ncol(counts)], na.rm = TRUE, dims = 1) # nolint

total_counts[1]

nm_fpkm <- t(do.call(
    rbind,
    lapply(
        1:length(total_counts), # nolint
        function(i) {
            10^9 * counts[, i] / lens_num / total_counts[i] # nolint
        }
    )
))
counts_colname <- colnames(counts)
colnames(nm_fpkm) <- counts_colname
write.csv(nm_fpkm, file = "GSE75246.csv")

# 文档结束————————————————————————————————————————————————————————————————————————————————————————————————————



counts <- read.delim("GSE79819merge.all_RNA_counts.txt",
    header = TRUE,
    fill = TRUE,
    row.names = 1
)
len <- read.delim("GSE79819merge.all_RNA_lengths.txt",
    header = TRUE,
    fill = TRUE
)
lens_name <- len$ID
lengths <- as.data.frame(lapply( # nolint
    len, # nolint
    as.numeric
))
lens_num <- select(lengths, 2)
lens_num <- unlist(lens_num)

# ng <- intersect(rownames(counts), len$ID)
# length(ng)
# lens <- len[match(ng, len$ID), 2]
names(lens_num) <- lens_name
head(lens_num)

counts[, 1:ncol(counts)] <- as.data.frame(lapply( # nolint
    counts[, 1:ncol(counts)], # nolint
    as.numeric
))

total_counts <- colSums(counts[, 1:ncol(counts)], na.rm = TRUE, dims = 1) # nolint

total_counts[1]

nm_fpkm <- t(do.call(
    rbind,
    lapply(
        1:length(total_counts), # nolint
        function(i) {
            10^9 * counts[, i] / lens_num / total_counts[i] # nolint
        }
    )
))
counts_colname <- colnames(counts)
colnames(nm_fpkm) <- counts_colname
write.csv(nm_fpkm, file = "GSE79819.csv")

# 文档结束————————————————————————————————————————————————————————————————————————————————————————————————————



counts <- read.delim("GSE99074merge.all_RNA_counts.txt",
    header = TRUE,
    fill = TRUE,
    row.names = 1
)
len <- read.delim("GSE99074merge.all_RNA_lengths.txt",
    header = TRUE,
    fill = TRUE
)
lens_name <- len$ID
lengths <- as.data.frame(lapply( # nolint
    len, # nolint
    as.numeric
))
lens_num <- select(lengths, 2)
lens_num <- unlist(lens_num)

# ng <- intersect(rownames(counts), len$ID)
# length(ng)
# lens <- len[match(ng, len$ID), 2]
names(lens_num) <- lens_name
head(lens_num)

counts[, 1:ncol(counts)] <- as.data.frame(lapply( # nolint
    counts[, 1:ncol(counts)], # nolint
    as.numeric
))

total_counts <- colSums(counts[, 1:ncol(counts)], na.rm = TRUE, dims = 1) # nolint

total_counts[1]

nm_fpkm <- t(do.call(
    rbind,
    lapply(
        1:length(total_counts), # nolint
        function(i) {
            10^9 * counts[, i] / lens_num / total_counts[i] # nolint
        }
    )
))
counts_colname <- colnames(counts)
colnames(nm_fpkm) <- counts_colname
write.csv(nm_fpkm, file = "GSE99074.csv")

# 文档结束————————————————————————————————————————————————————————————————————————————————————————————————————