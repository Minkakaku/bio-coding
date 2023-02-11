library(Seurat)
library(R.utils)
library(dplyr)
library(Matrix)
list.files()
# load data (data & annote)
data <- readRDS("data/GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds")
annote <- read.delim(
    file = "data/GSE131907_Lung_Cancer_cell_annotation.txt",
    header = TRUE
)
annote <- annote[, c(7:7)]
# process data
# data <- dplyr::filter(data, !duplicated(data$X))
# rownames(data) <- data[, 1]
# data <- data[, -1]
cellname <- paste(rep("C", times = ncol(data)),
    seq_len(ncol(data)),
    sep = "_"
)
colnames(data) <- cellname
rownames(annote) <- colnames(data)
annote <- cbind(allele = row.names(annote), annote)
row.names(annote) <- NULL

colnames(annote) <- c("", "Cell", "Cell_type")
# load NCBI gene information
geneinfo <- read.delim("/data/gpfs02/wmo/work/HF/scRNA-seq/HF/220731scCATCH/Homo_sapiens.gene_info")
geneinfo <- geneinfo[, c(3, 5)]
geneinfo$species <- rep("Human", times = nrow(geneinfo))
# create SeuratObject
data <- CreateSeuratObject(counts = data)
# data <- NormalizeData(object = data)
# generate data.rds
# saveRDS(data, file = "data/data.rds")

# data is Seurat object after log-normalization
human_trainsets <- readRDS("data/data.rds")
human_trainsets <- human_trainsets[["RNA"]]@data

# revising gene symbols
genename <- rownames(human_trainsets)
genename1 <- genename[genename %in% geneinfo$Symbol]
genename2 <- genename[!genename %in% geneinfo$Symbol]
genename3 <- genename2[genename2 %in% geneinfo$Synonyms]
genename4 <- rep("NA", length(genename3))
for (i in seq_len(length(genename3))) {
    d1 <- geneinfo[geneinfo$Synonyms == genename3[i], ]$Symbol
    if (length(d1) == 1) {
        genename4[i] <- d1
    }
}
genename3 <- c(genename1, genename3)
genename4 <- c(genename1, genename4)
genedata <- data.frame(raw_name = genename3, new_name = genename4, stringsAsFactors = F)
genedata <- genedata[!genedata$new_name == "NA", ]
genedata1 <- as.data.frame(table(genedata$new_name), stringsAsFactors = F)
genedata1 <- genedata1[genedata1$Freq == 1, ]
genedata <- genedata[genedata$new_name %in% genedata1$Var1, ]
human_trainsets <- human_trainsets[genedata$raw_name, ]
all(rownames(human_trainsets) == genedata$raw_name)
rownames(human_trainsets) <- genedata$new_name
all(rownames(human_trainsets) == genedata$new_name)
all(rownames(human_trainsets) %in% geneinfo$Symbol)

# for log-normalization, please execute the following command
# human_trainsets<- CreateSeuratObject(counts = human_trainsets)
# human_trainsets<- NormalizeData(object = human_trainsets)
# human_trainsets<- human_trainsets[['RNA']]@data

library(Rcpp)
Rcpp::sourceCpp(code = "
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix asMatrix(NumericVector rp,
                       NumericVector cp,
                       NumericVector z,
                       int nrows,
                       int ncols){
  int k = z.size() ;
  NumericMatrix  mat(nrows, ncols);
  for (int i = 0; i < k; i++){
      mat(rp[i],cp[i]) = z[i];
  }
  return mat;
}
")

as.matrix <- function(mat) {
    row_pos <- mat@i
    col_pos <- findInterval(seq(mat@x) - 1, mat@p[-1])
    tmp <- asMatrix(
        rp = row_pos, cp = col_pos, z = mat@x,
        nrows = mat@Dim[1], ncols = mat@Dim[2]
    )
    row.names(tmp) <- mat@Dimnames[[1]]
    colnames(tmp) <- mat@Dimnames[[2]]
    return(tmp)
}

human_trainsets <- as.matrix(human_trainsets)
write.csv(human_trainsets, file = "data/human_trainsets.csv")
gzip("data/human_trainsets.csv", "data/human_trainsets.csv.gz")
write.csv(annote, file = "data/annote.csv")

# human_trainsets_data.csv or human_trainsets_data.gz can be used for running scDeepSort

df = read.csv("/data/gpfs02/wmo/work/HF/scRNA-seq/HF/220712DeepClassfier/data/annote.csv")
annote <- df[, c(2, 3)]
head(annote)
colnames(annote)
