rm(list = ls())
# hongfan
library(KEGGREST)

# return related genes
rege <- keggFind("genes", "complement")
rege2 <- keggFind("genes", "DNA binding protein")
rege3 <- keggFind("genes", "Oligoadenylate synthetase")

# Find genes
gene <- character()
for (i in seq_len(length(rege))) {
    gene1 <- unlist(strsplit(rege[i], ";"))[1]
    if (unlist(strsplit(names(gene1), ":"))[1] == "mmu") {
        gene <- cbind(gene, gene1)
    }
}
genes <- unlist(lapply(gene, function(x) strsplit(x, ",")))
library(stringr)
genes <- str_trim(genes, "both")
genes1 <- unique(genes)

dir <- "/home/hongfan/XMU2021/JOBs/gene_process/yhf/4.Differential/1.deglist/group2vsgroup1/"
new <- read.csv(paste0(dir, "group2vsgroup1_deg.csv"))
new <- new[new$gene_name %in% genes1, ]
write.csv(new, "/home/hongfan/XMU2021/JOBs/gene_process/yhf/complement related genes.csv")
