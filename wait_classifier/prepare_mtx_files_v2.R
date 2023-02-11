rm(list = ls())

library(dplyr)
library(Matrix)
library(Seurat)


assays <- dir('/home/wsa/data/scRNAlist_1850_6590_h5seurat')
dir <- paste0("/home/wsa/data/scRNAlist_1850_6590_h5seurat/", assays)
seurat.1 <- SeuratDisk::LoadH5Seurat(dir[[1]])
seurat.2 <- SeuratDisk::LoadH5Seurat(dir[[2]])
seurat.3 <- SeuratDisk::LoadH5Seurat(dir[[3]])
seurat.4 <- SeuratDisk::LoadH5Seurat(dir[[4]])
seurat.5 <- SeuratDisk::LoadH5Seurat(dir[[5]])
seurat.6 <- SeuratDisk::LoadH5Seurat(dir[[6]])
seurat.7 <- SeuratDisk::LoadH5Seurat(dir[[7]])


seurat <- merge(x=seurat.1, y=list(seurat.2,seurat.3,seurat.5,seurat.7,seurat.4, seurat.6 ),all = FALSE)
#seurat <- subset(seurat, cells=WhichCells(seurat, expression = Tissue %in% c('Blood','Lung','Liver','Spleen')))
ncol(seurat)
table(seurat@meta.data$group)
table(seurat@meta.data$cluster)
seurat@meta.data$cluster[which(seurat@meta.data$cluster =='TAM 1')] <- 'MDM'
seurat@meta.data$cluster[which(seurat@meta.data$cluster =='TAM 2')] <- 'Microglia'
seurat@meta.data$cluster[which(seurat@meta.data$cluster =='DC')] <- 'Dendritic Cells'
seurat@meta.data$cluster[which(seurat@meta.data$cluster =='prol. TAM')] <- 'TAM 3'
#table(seurat@meta.data$cluster)

seurat <- NormalizeData(seurat)

#seurat$cluster <- plyr::mapvalues(seurat$cluster,
#                                  from = c('G3','G4','G5a','G5b','G5c'),
#                                  to = c('immNeu','mNeu','PMNa','PMNb','PMNc'))


################################################################################
## 对每一个课题组的数据找高变基因，然后再合起来

seurat.list <- list(seurat.1, 
                    seurat.2,
                    seurat.3,
                    seurat.4, 
                    seurat.5,
                    seurat.6,
                    seurat.7
                    )

hvgs.list <- list()
cnt <- 1
for (ss in seurat.list) {
  ss <- NormalizeData(ss)
  ss <- FindVariableFeatures(ss, nfeatures=2000)
  hvgs.list[[cnt]] <- VariableFeatures(ss)
  cnt <- cnt+1
}

hvgs <- hvgs.list %>% purrr::reduce(union)
#hvgs <- union(hvgs, c('Ltb4r1','Cxcr2','Cxcr4','Cd274'))
cat(length(hvgs))


genes.list <- lapply(seurat.list, rownames)
genes <- genes.list %>% purrr::reduce(intersect)
hvgs <- intersect(hvgs, genes)

###############################################################################
##将有标签anno的数据保存起来,mtx
# cells <- colnames(seurat.xie2020.ctrl)
# cells <- WhichCells(seurat.xie2020.ctrl, expression= Tissue!='Bone marrow')
cells <- WhichCells(seurat, expression= (group %in% c(4,6)))
length(cells)
X <- GetAssayData(seurat, 'data')
X <- t(X[hvgs, ] %>% as.matrix)
X[X<0] <- 0

cat('write data\n')
setwd('/home/wsa/data/mtx_file_1850_6590')
spMat <- Matrix(X[cells, ], sparse = TRUE)
cat(nrow(spMat), ',', ncol(spMat), '\n')
writeMM(spMat, file = paste('ref_dataset.mtx',sep=''))

# save label (text)
write.table(seurat@meta.data[cells, 'cluster'] %>% as.character,
            file = paste('ref_dataset_cluster.txt',sep=''),
            sep = '\n', row.names = F, col.names = F, quote = T)

write.table(hvgs, file = 'ref_dataset_genes.txt', row.names = F, col.names = F)
write.table(cells, file = 'ref_dataset_cells.txt', row.names = F, col.names = F)

################################################################################
## 将没有标签anno的数据保存起来,mtx
# seurat.alt <- merge(x=seurat.ni2022, y=list(seurat.xie2020.eco, seurat.immgen))
# cells <- colnames(seurat.alt)
cells0 <- WhichCells(seurat)
cells <- setdiff(cells0, cells)
length(cells)

X <- GetAssayData(seurat, 'data')
X <- t(X[hvgs, ] %>% as.matrix)
X[X<0] <- 0

cat('write data\n')
spMat <- Matrix(X[cells, ], sparse = TRUE)
cat(nrow(spMat), ',', ncol(spMat), '\n')
writeMM(spMat, file = paste('alt_dataset.mtx',sep=''))


write.table(hvgs, file = 'alt_dataset_genes.txt', row.names = F, col.names = F)
write.table(cells, file = 'alt_dataset_cells.txt', row.names = F, col.names = F)



################################################################################
## 将所有数据保存起来,mtx
# seurat <- merge(x=seurat.ni2022, y=list(seurat.xie2020.ctrl, seurat.xie2020.eco, seurat.immgen))
cells <- WhichCells(seurat)


X <- GetAssayData(seurat, 'data')
X <- t(X[hvgs, ] %>% as.matrix)
X[X<0] <- 0

cat('write data\n')
spMat <- Matrix(X[cells, ], sparse = TRUE)
cat(nrow(spMat), ',', ncol(spMat), '\n')
writeMM(spMat, file = paste('all_dataset.mtx',sep=''))


write.table(hvgs, file = 'all_dataset_genes.txt', row.names = F, col.names = F)
write.table(cells, file = 'all_dataset_cells.txt', row.names = F, col.names = F)



################################################################################
## 保存合并后的seurat对象为h5seurat格式
setwd('/home/wsa/data/')
seurat <- seurat[,cells]

SeuratDisk::SaveH5Seurat(seurat, filename = 'all_1850_6590.h5seurat', overwrite = TRUE)
#SeuratDisk::Convert('all3.h5seurat', dest = 'all_1850_6590.h5ad', overwrite = TRUE)

################################################################################
## 保存用于伪时序分析的seurat对象
alt_cell_type <- read.table('/home/wsa/data/mtx_file_1850_6590/alt_cluster_predict.txt',sep = ',', header = T, row.names = 1,col.names = "V1")
ref_cell_type <- read.table('/home/wsa/data/mtx_file_1850_6590/ref_dataset_cluster.txt',sep = '\n')
cell_type = rbind(alt_cell_type,ref_cell_type)
cell <- Cells(seurat)
cells_type = cbind(cell,cell_type)
rownames(cells_type)=cells_type[,1]  #取出第一列

# AddMetaData一定要保证cell_type的行名和sdata的行名一样，否则会出现NA
seurat <- AddMetaData(
  object = seurat,
  metadata = cells_type,
  col.name = names(cells_type))

setwd('/home/wsa/data/')
cells <- WhichCells(seurat, expression= (V1 %in% c('TAM 1','TAM 2','Monocytes','prol. TAM')))
seurat <- seurat[,cells]
SeuratDisk::SaveH5Seurat(seurat, filename = 'TI_1850_6590.h5seurat', overwrite = TRUE)
