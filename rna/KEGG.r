rm(list = ls())
# hongfan
library(KEGGREST)

# return related genes
rege <- keggFind("genes", "RNA binding protein")
rege2 <- keggFind("genes", "DNA binding protein")
rege3 <- keggFind("genes", "Oligoadenylate synthetase")

# Find genes
gene <- character()
for (i in seq_len(length(rege3))) {
    gene1 <- unlist(strsplit(rege3[i], ";"))[1]
    if (unlist(strsplit(names(gene1), ":"))[1] == "mmu") {
        gene <- cbind(gene, gene1)
    }
}
genes <- unlist(lapply(gene, function(x) strsplit(x, ",")))
genes1 <- unique(genes)
