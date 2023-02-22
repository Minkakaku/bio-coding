setwd(".")
path <- "."

multmerge <- function(path) {
    filenames <- list.files(
        path = path,
        pattern = ".counts.txt"
    )
    datalist <- lapply(filenames, function(x) {
        read.delim(
            file = x, header = TRUE,
            col.names = c("geneid", unlist(strsplit(x, ".", fixed = TRUE))[1])
        )
    })
    Reduce(function(x, y) {
        merge(x, y,
            by = "geneid"
        )
    }, datalist)
}

mergedata <- multmerge(path)
write.csv(mergedata, "merge.counts.csv")
