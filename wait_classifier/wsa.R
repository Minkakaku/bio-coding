rm(list = ls())

library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratDisk)

scRNAlist_tmp <- list()
scRNAlist <- list()
################################################################################
#1111111111读取GSE117891，一个count.tsv 和一个注释文件
################################################################################
#读取矩阵
setwd('/home/wsa/data/young_older/GSE117891')
assays <- dir('/home/wsa/data/young_older/GSE117891')
dir <- paste0("/home/wsa/data/young_older/GSE117891/", assays)
counts <- read.table(dir, header=TRUE, row.names=1)
ncol(counts)

#读取临床信息
annotation = read.table("/home/wsa/data/young_older_clinical/GSE117891/GSE117891_Sample_barcode_cell_information.txt",header = T,sep='\t',row.names = 1)
nrow(annotation)
clinical = read.table("/home/wsa/data/young_older_clinical/GSE117891/clinical.txt",header = T,sep = "\t",fileEncoding = 'utf-8')
col = colnames(clinical)
sample = clinical$sample
for (j in sample){
  anno_rownumber <- which(annotation$sample == j, arr.ind = TRUE)
  cli_rownumber <- which(clinical$sample == j, arr.ind = TRUE)
  annotation[anno_rownumber,paste(col)] <- clinical[cli_rownumber,]
}
nrow(annotation)

#添加年龄分级信息
anno_rownumber <- which(annotation$age <= 50 &annotation$age >=18, arr.ind = TRUE)
annotation[anno_rownumber,'level'] <- 'adult'
anno_rownumber <- which(annotation$age <= 100 & annotation$age >=65, arr.ind = TRUE)
annotation[anno_rownumber,'level'] <- 'aged'
View(annotation)

#将矩阵和metadata变成seuratobject
scRNAlist_tmp[[1]] <- CreateSeuratObject(counts,min.cells=3, min.features = 200)
scRNAlist_tmp[[1]] <- AddMetaData(
  object = scRNAlist_tmp[[1]],
  metadata = annotation,
  col.name = c(names(annotation),"note"))

#scRNAlist_tmp[[1]]是所有数据，下面选出我想要的数据
tmp <- subset(scRNAlist_tmp[[1]]@meta.data, tumor %in% c("GBM") & level %in% c("adult","aged"))
scRNAlist[[1]] <- subset(scRNAlist_tmp[[1]], cells=row.names(tmp))
scRNAlist[[1]] <- RenameCells(scRNAlist[[1]], add.cell.id = 1)
rm(tmp)

if(FALSE){
################################################################################
#读取GSE173278，三个标准文件加一个metadata文件
################################################################################
#读取矩阵
setwd('/home/wsa/data/young_older/GSE173278')
assays <- dir('/home/wsa/data/young_older/GSE173278')
dir <- paste0("/home/wsa/data/young_older/GSE173278/", assays)
gene_file = read.table(dir[[2]], stringsAsFactors=F, header=F)
if (dim(gene_file)[2] == 1){
  genes = read.table(dir[[2]], stringsAsFactors=F, header=F)$V1
  genes = make.unique(genes, sep = '.')
} else
{
  genes = read.table(dir[[2]], stringsAsFactors=F, header=F)$V2
  genes = make.unique(genes, sep = '.')
}
barcodes = readLines(dir[[1]])
mtx = Matrix::readMM(dir[[3]])
mtx = as(mtx, 'dgCMatrix')
colnames(mtx) = barcodes
rownames(mtx) = genes
ncol(mtx)

#读取临床信息
annotation = read.table("/home/wsa/data/young_older_clinical/GSE173278/meta.csv",header = T,sep=',',row.names = 1)
clinical = read.table("/home/wsa/data/young_older_clinical/GSE173278/clinical2.txt",header = T,sep = "\t",fileEncoding = 'utf-8')
col = colnames(clinical)
sample = clinical$sample
for (j in sample){
  anno_rownumber <- which(annotation$sample == j, arr.ind = TRUE)
  cli_rownumber <- which(clinical$sample == j, arr.ind = TRUE)
  annotation[anno_rownumber,paste(col)] <- clinical[cli_rownumber,]
}
nrow(annotation)
#添加年龄分级信息
anno_rownumber <- which(annotation$age <= 50, arr.ind = TRUE)
annotation[anno_rownumber,'level'] <- 'adult'
anno_rownumber <- which(annotation$age <= 70 & annotation$age >=60, arr.ind = TRUE)
annotation[anno_rownumber,'level'] <- 'aged'

#将矩阵和metadata变成seuratobject
scRNAlist_tmp[[2]] <- CreateSeuratObject(mtx,min.cells=3, min.features = 200)
scRNAlist_tmp[[2]] <- AddMetaData(
  object = scRNAlist_tmp[[2]],
  metadata = annotation,
  col.name = names(annotation))
ncol(scRNAlist_tmp[[2]])
#scRNAlist_tmp[[2]]是所有数据，下面选出我想要的数据
tmp <- subset(scRNAlist_tmp[[2]]@meta.data, sample %in% c("JK124","JK142","JK152","JK153","JK156","JK163"))
tis_org <- c(unique(tmp$note))
tis <- c( "JK124_reg1_tis_1",      "JK124_reg1_tis_2" ,  "JK124_reg2_tis_1"  ,    "JK124_reg2_tis_2" ,    "JK142_reg1_tis_1",   
          "JK142_reg1_cell_1",     "JK142_reg2_tis_1",      "JK142_reg2_tis_2.1_br", "JK142_reg2_tis_2.2_br",
         "JK142_reg2_cell_1",     
          "JK152_reg1_tis_1",     
            "JK152_reg1_cell_1",     "JK152_reg2_tis_1",           
                "JK152_reg2_cell_1",     "JK153_reg1_tis_1",        
          "JK153_reg1_cell_1",     "JK153_reg2_tis_1",          
          "JK156_reg1_tis_1",            "JK156_reg2_tis_1",     
          "JK156_reg2_tis_2_br",        "JK163_reg1_tis_1",           
              "JK163_reg1_cell_1",     "JK163_reg2_tis_1",           
          "JK163_reg2_cell_1")
tmp <- subset(tmp, note %in% tis)
nrow(tmp)
scRNAlist[[2]] <- subset(scRNAlist_tmp[[2]], cells=row.names(tmp))
scRNAlist[[2]] <- RenameCells(scRNAlist[[2]], add.cell.id = 2)
ncol(scRNAlist[[2]])
}
################################################################################
#222222读取GSE131928，一个tsv文件，10x
################################################################################
#读取矩阵
setwd('/home/wsa/data/young_older/GSE131928')
assays <- dir('/home/wsa/data/young_older/GSE131928')
dir <- paste0("/home/wsa/data/young_older/GSE131928/", assays)
counts2 <- read.table(dir[[2]], header=TRUE, row.names=1,check.names=FALSE)
ncol(counts2)

#读取临床信息
annotation = read.table("/home/wsa/data/young_older_clinical/GSE131928/anno8.txt",header = T,sep='\t',row.names = 1)
clinical = read.table("/home/wsa/data/young_older_clinical/GSE131928/clinical8.txt",header = T,sep = "\t",fileEncoding = 'utf-8')
col = colnames(clinical)
sample = clinical$sample
for (j in sample){
  anno_rownumber <- which(annotation$sample %in% j, arr.ind = TRUE)
  cli_rownumber <- which(clinical$sample %in% j, arr.ind = TRUE)
  annotation[anno_rownumber,paste(col)] <- clinical[cli_rownumber,]
}

#添加年龄分级信息
anno_rownumber <- which(annotation$age <= 50 &annotation$age >=18, arr.ind = TRUE)
annotation[anno_rownumber,'level'] <- 'adult'
anno_rownumber <- which(annotation$age <= 100 & annotation$age >=65, arr.ind = TRUE)
annotation[anno_rownumber,'level'] <- 'aged'
View(annotation)

#将矩阵和metadata变成seuratobject
scRNAlist_tmp[[2]] <- CreateSeuratObject(counts2,min.cells=3, min.features = 200)
scRNAlist_tmp[[2]] <- AddMetaData(
  object = scRNAlist_tmp[[2]],
  metadata = annotation,
  col.name = names(annotation))

#scRNAlist_tmp[[3]]是所有数据，下面选出我想要的数据
tmp <- subset(scRNAlist_tmp[[2]]@meta.data, tumor %in% c("GBM") & level %in% c("adult","aged") & new %in% c("Primary","primary"))
scRNAlist[[2]] <- subset(scRNAlist_tmp[[2]], cells=row.names(tmp))
scRNAlist[[2]] <- RenameCells(scRNAlist[[2]], add.cell.id = 2)
rm(tmp)


################################################################################
## 33333333读取TCGA，三个标准文件
################################################################################
assays <- dir("/home/wsa/data/young_older/TCGA/")
dir <- paste0("/home/wsa/data/young_older/TCGA/", assays)
cli_data = read.table("/home/wsa/data/young_older_clinical/TCGA/clinical.txt",header = T,sep = "\t")
scRNAtmp <- list()
for(i in seq(1,length(dir),3)){
  genes = read.table(dir[[i+1]], stringsAsFactors=F, header=F)$V2
  genes = make.unique(genes, sep = '.')
  barcodes = readLines(dir[[i]])
  mtx = Matrix::readMM(dir[[i+2]])
  mtx = as(mtx, 'dgCMatrix')
  colnames(mtx) = barcodes
  rownames(mtx) = genes
  
  #不设置min.cells过滤基因会导致CellCycleScoring报错：
  #Insufficient data values to produce 24 bins.  
  scRNAtmp[[i%/%3+1]] <- CreateSeuratObject(mtx ,min.cells=3, min.features = 200)
  
  #添加临床信息
  scRNAtmp[[i%/%3+1]] <- AddMetaData(
    object = scRNAtmp[[i%/%3+1]],
    metadata = c(rep(cli_data[i%/%3+1,],scRNAtmp[[i%/%3+1]]@assays$RNA@counts@Dim[2])),
    col.name = colnames(cli_data)
  )
  #每一组样本的细胞都加上标签
  scRNAtmp[[i%/%3+1]] <- RenameCells(scRNAtmp[[i%/%3+1]], 
                                     new.names = paste0(Cells(x = scRNAtmp[[i%/%3+1]]), i%/%3+1 ))
}

scRNAlist[[3]] <- merge(x=scRNAtmp[[1]], y=list(scRNAtmp[[2]]),all = FALSE)
#tmp <- subset(scRNAlist[[8]]@meta.data, level %in% c("adult_normal","aged_normal"))
#scRNAlist[[8]] <- subset(scRNAlist[[8]], cells=row.names(tmp))
#给细胞barcode加个前缀，防止合并后barcode重名
scRNAlist[[3]] <- RenameCells(scRNAlist[[3]], add.cell.id = 3)


if (FALSE) {
################################################################################
#333333333读取GSE131928，一个tsv文件，smartseq
################################################################################
#读取矩阵
setwd('/home/wsa/data/young_older/GSE131928')
assays <- dir('/home/wsa/data/young_older/GSE131928')
dir <- paste0("/home/wsa/data/young_older/GSE131928/", assays)
counts1 <- read.table(dir[[1]], header=TRUE, row.names=1,check.names=FALSE)
ncol(counts1)
    
#读取临床信息
annotation = read.table("/home/wsa/data/young_older_clinical/GSE131928/anno3.txt",header = T,sep='\t',row.names = 1)
clinical = read.table("/home/wsa/data/young_older_clinical/GSE131928/clinical3.txt",header = T,sep = "\t",fileEncoding = 'utf-8')
col = colnames(clinical)
sample = clinical$sample
for (j in sample){
  anno_rownumber <- which(annotation$sample == j, arr.ind = TRUE)
  cli_rownumber <- which(clinical$sample == j, arr.ind = TRUE)
  annotation[anno_rownumber,paste(col)] <- clinical[cli_rownumber,]
}

#添加年龄分级信息
anno_rownumber <- which(annotation$age <= 50 &annotation$age >=18, arr.ind = TRUE)
annotation[anno_rownumber,'level'] <- 'adult'
anno_rownumber <- which(annotation$age <= 100 & annotation$age >=65, arr.ind = TRUE)
annotation[anno_rownumber,'level'] <- 'aged'
View(annotation)

#将矩阵和metadata变成seuratobject
scRNAlist_tmp[[3]] <- CreateSeuratObject(counts1,min.cells=3, min.features = 200)
scRNAlist_tmp[[3]] <- AddMetaData(
  object = scRNAlist_tmp[[3]],
  metadata = annotation,
  col.name = names(annotation))

#scRNAlist_tmp[[3]]是所有数据，下面选出我想要的数据
tmp <- subset(scRNAlist_tmp[[3]]@meta.data, tumor %in% c("GBM") & level %in% c("adult","aged") & new %in% c("Primary","primary"))
scRNAlist[[3]] <- subset(scRNAlist_tmp[[3]], cells=row.names(tmp))
scRNAlist[[3]] <- RenameCells(scRNAlist[[3]], add.cell.id = 3)
rm(tmp)
ncol(scRNAlist[[3]])
}

################################################################################
#4444444444读取GSE163120，一个csv文件
################################################################################
#读取矩阵
setwd('/home/wsa/data/young_older/GSE163120')
assays <- dir('/home/wsa/data/young_older/GSE163120')
dir <- paste0("/home/wsa/data/young_older/GSE163120/", assays)
counts <- read.table(dir, header=TRUE, row.names=1,sep = ',',check.names = FALSE)
ncol(counts)

#读取临床信息
annotation = read.table("/home/wsa/data/young_older_clinical/GSE163120/anno4.txt",header = T,sep='\t',row.names = 1)
clinical = read.table("/home/wsa/data/young_older_clinical/GSE163120/clinical4.txt",header = T,sep = "\t",fileEncoding = 'utf-8')
col = colnames(clinical)
sample = clinical$sample
for (j in sample){
  anno_rownumber <- which(annotation$sample == j, arr.ind = TRUE)
  cli_rownumber <- which(clinical$sample == j, arr.ind = TRUE)
  annotation[anno_rownumber,paste(col)] <- clinical[cli_rownumber,]
}

#添加年龄分级信息
anno_rownumber <- which(annotation$age <= 50 &annotation$age >=18, arr.ind = TRUE)
annotation[anno_rownumber,'level'] <- 'adult'
anno_rownumber <- which(annotation$age <= 100 & annotation$age >=65, arr.ind = TRUE)
annotation[anno_rownumber,'level'] <- 'aged'
View(annotation)

#将矩阵和metadata变成seuratobject
scRNAlist_tmp[[4]] <- CreateSeuratObject(counts,min.cells=3, min.features = 200)
scRNAlist_tmp[[4]] <- AddMetaData(
  object = scRNAlist_tmp[[4]],
  metadata = annotation,
  col.name = names(annotation))

#scRNAlist_tmp[[1]]是所有数据，下面选出我想要的数据
tmp <- subset(scRNAlist_tmp[[4]]@meta.data, tumor %in% c("GBM") & level %in% c("adult","aged") & new %in% c("Primary","primary"))
scRNAlist[[4]] <- subset(scRNAlist_tmp[[4]], cells=row.names(tmp))
scRNAlist[[4]] <- RenameCells(scRNAlist[[4]], add.cell.id = 4)
ncol(scRNAlist[[4]])
rm(tmp)

if(FALSE){
################################################################################
#55555555读取GSE102130，一个csv文件,无annotation
################################################################################
#读取矩阵
setwd('/home/wsa/data/young_older/GSE102130')
assays <- dir('/home/wsa/data/young_older/GSE102130')
dir <- paste0("/home/wsa/data/young_older/GSE102130/", assays)
counts <- read.table(dir, header=TRUE, row.names=1,sep = '\t',check.names = FALSE)
ncol(counts)

#读取临床信息
annotation <- list()
topic <- strsplit(colnames(counts), "-", fixed= T)
for (i in 1:length(topic)){
  annotation[[i]] = topic[[i]][1]
}
annotation <- data.frame(matrix(unlist(annotation), nrow=length(topic)))
names(annotation) = c('sample')
row.names(annotation) = colnames(counts)
clinical = read.table("/home/wsa/data/young_older_clinical/GSE102130/clinical.txt",header = T,sep = "\t",fileEncoding = 'utf-8')
col = colnames(clinical)
sample = clinical$sample
for (j in sample){
  anno_rownumber <- which(annotation$sample == j, arr.ind = TRUE)
  cli_rownumber <- which(clinical$sample == j, arr.ind = TRUE)
  annotation[anno_rownumber,paste(col)] <- clinical[cli_rownumber,]
}

#添加年龄分级信息
anno_rownumber <- which(annotation$age <= 50 &annotation$age >=18, arr.ind = TRUE)
annotation[anno_rownumber,'level'] <- 'adult'
anno_rownumber <- which(annotation$age <= 100 & annotation$age >=65, arr.ind = TRUE)
annotation[anno_rownumber,'level'] <- 'aged'
View(annotation)

#将矩阵和metadata变成seuratobject
scRNAlist_tmp[[5]] <- CreateSeuratObject(counts,min.cells=3, min.features = 200)
scRNAlist_tmp[[5]] <- AddMetaData(
  object = scRNAlist_tmp[[5]],
  metadata = annotation,
  col.name = names(annotation))

#scRNAlist_tmp[[1]]是所有数据，下面选出我想要的数据
tmp1 <- subset(scRNAlist_tmp[[5]]@meta.data, age<=50)
tmp2 <- subset(scRNAlist_tmp[[5]]@meta.data, age<=70)
tmp2 <- subset(tmp2, age>=60)
scRNAlist[[5]] <- subset(scRNAlist_tmp[[5]], cells=c(row.names(tmp1),row.names(tmp2)))
scRNAlist[[5]] <- RenameCells(scRNAlist[[5]], add.cell.id = 5)
ncol(scRNAlist[[5]])
}

################################################################################
#55555555读取GSE141946，一个txt文件
################################################################################
#读取矩阵
setwd('/home/wsa/data/young_older/GSE141946')
assays <- dir('/home/wsa/data/young_older/GSE141946')
dir <- paste0("/home/wsa/data/young_older/GSE141946/", assays)
counts <- read.table(dir, header=TRUE, row.names=1,sep = '\t',check.names = FALSE)
ncol(counts)

#读取临床信息
annotation = read.table("/home/wsa/data/young_older_clinical/GSE141946/anno_6.txt",header = T,sep='\t',row.names = 1)
clinical = read.table("/home/wsa/data/young_older_clinical/GSE141946/clinical_6.txt",header = T,sep = "\t",fileEncoding = 'utf-8')
col = colnames(clinical)
sample = clinical$sample
for (j in sample){
  anno_rownumber <- which(annotation$sample == j, arr.ind = TRUE)
  cli_rownumber <- which(clinical$sample == j, arr.ind = TRUE)
  annotation[anno_rownumber,paste(col)] <- clinical[cli_rownumber,]
}


#添加年龄分级信息
anno_rownumber <- which(annotation$age <= 50 &annotation$age >=18, arr.ind = TRUE)
annotation[anno_rownumber,'level'] <- 'adult'
anno_rownumber <- which(annotation$age <= 100 & annotation$age >=65, arr.ind = TRUE)
annotation[anno_rownumber,'level'] <- 'aged'
View(annotation)

#将矩阵和metadata变成seuratobject
scRNAlist_tmp[[5]] <- CreateSeuratObject(counts,min.cells=3, min.features = 200)
scRNAlist_tmp[[5]] <- AddMetaData(
  object = scRNAlist_tmp[[5]],
  metadata = annotation,
  col.name = names(annotation))

#scRNAlist_tmp[[1]]是所有数据，下面选出我想要的数据
tmp <- subset(scRNAlist_tmp[[5]]@meta.data, tumor %in% c("GBM") & level %in% c("adult","aged") & new %in% c("Primary","primary"))
scRNAlist[[5]] <- subset(scRNAlist_tmp[[5]], cells=row.names(tmp))
scRNAlist[[5]] <- RenameCells(scRNAlist[[5]], add.cell.id = 5)
ncol(scRNAlist[[5]])
rm(tmp)

if(FALSE) {
  
################################################################################
#读取GSE70630，一个txt文件,无annotation
################################################################################
#读取矩阵
setwd('/home/wsa/data/young_older/GSE70630')
assays <- dir('/home/wsa/data/young_older/GSE70630')
dir <- paste0("/home/wsa/data/young_older/GSE70630/", assays)
counts <- read.table(dir, header=TRUE, row.names=1,sep = '\t',check.names = FALSE)

#读取临床信息
annotation <- list()
topic <- strsplit(colnames(counts), "_", fixed= T)
for (i in 1:length(topic)){
  annotation[[i]] = topic[[i]][1]
}
annotation <- data.frame(matrix(unlist(annotation), nrow=length(topic)))
names(annotation) = c('sample')
row.names(annotation) = colnames(counts)
clinical = read.table("/home/wsa/data/young_older_clinical/GSE70630/clinical7.txt",header = T,sep = "\t",fileEncoding = 'utf-8')
col = colnames(clinical)
sample = clinical$sample
for (j in sample){
  anno_rownumber <- which(annotation$sample == j, arr.ind = TRUE)
  cli_rownumber <- which(clinical$sample == j, arr.ind = TRUE)
  annotation[anno_rownumber,paste(col)] <- clinical[cli_rownumber,]
}

#添加年龄分级信息
anno_rownumber <- which(annotation$age <= 50, arr.ind = TRUE)
annotation[anno_rownumber,'level'] <- 'young'
anno_rownumber <- which(annotation$age <= 70 & annotation$age >=60, arr.ind = TRUE)
annotation[anno_rownumber,'level'] <- 'older'

#将矩阵和metadata变成seuratobject
scRNAlist_tmp[[6]] <- CreateSeuratObject(counts,min.cells=3, min.features = 200)
scRNAlist_tmp[[6]] <- AddMetaData(
  object = scRNAlist_tmp[[6]],
  metadata = annotation,
  col.name = names(annotation))

#scRNAlist_tmp[[1]]是所有数据，下面选出我想要的数据
tmp1 <- subset(scRNAlist_tmp[[6]]@meta.data, age<=50)
tmp2 <- subset(scRNAlist_tmp[[6]]@meta.data, age<=70)
tmp2 <- subset(tmp2, age>=60)
scRNAlist[[6]] <- subset(scRNAlist_tmp[[6]], cells=c(row.names(tmp1),row.names(tmp2)))
scRNAlist[[6]] <- RenameCells(scRNAlist[[6]], add.cell.id = 6)

}

################################################################################
#6666666666读取group6，一个csv文件
################################################################################
#读取矩阵
setwd('/home/wsa/data/young_older/group6')
assays <- dir('/home/wsa/data/young_older/group6')
dir <- paste0("/home/wsa/data/young_older/group6/", assays)
counts <- read.table(dir, header=TRUE, row.names=1,sep = ',',check.names = FALSE)
ncol(counts)

#读取临床信息
annotation = read.table("/home/wsa/data/young_older_clinical/group6/anno_clinical7.txt",header = T,sep='\t',row.names = 1,)

#将矩阵和metadata变成seuratobject
scRNAlist_tmp[[6]] <- CreateSeuratObject(counts,min.cells=3, min.features = 200)
scRNAlist_tmp[[6]] <- AddMetaData(
  object = scRNAlist_tmp[[6]],
  metadata = annotation,
  col.name = names(annotation))

#加上名字前缀
scRNAlist[[6]] <-scRNAlist_tmp[[6]]
scRNAlist[[6]] <- RenameCells(scRNAlist[[6]], add.cell.id = 6)
scRNAlist[[6]] <- subset(scRNAlist[[6]], cells=WhichCells(scRNAlist[[6]], expression = cluster %in% c('Tumor','NormalBrain')))


################################################################################
#77777777读取GSE89567，一个txt文件
################################################################################
#读取矩阵
setwd('/home/wsa/data/young_older/group7')
assays <- dir('/home/wsa/data/young_older/group7')
dir <- paste0("/home/wsa/data/young_older/group7/", assays)
counts <- read.table(dir,header = TRUE, row.names=1,sep = '\t',check.names = FALSE)
View(counts)
ncol(counts)

#读取临床信息
annotation_tmp <- list()
topic <- strsplit(colnames(counts), "_", fixed= T)
for (i in 1:length(topic)){
  annotation_tmp[[i]] = topic[[i]][1]
}
annotation_tmp <- data.frame(matrix(unlist(annotation_tmp), nrow=length(topic)))
names(annotation_tmp) = c('sample')
annotation <- list()
topic2 <- strsplit(annotation_tmp$sample, "-", fixed= T)
for (i in 1:length(topic2)){
  annotation[[i]] = topic2[[i]][1]
}
annotation <- data.frame(matrix(unlist(annotation), nrow=length(topic)))
names(annotation) = c('sample')
row.names(annotation) = colnames(counts)
annotation[,c(1)][annotation[,c(1)] == "57"] = "MGH57"
annotation[,c(1)][annotation[,c(1)] == "MGH107neg"] = "MGH107"
annotation[,c(1)][annotation[,c(1)] == "MGH107pos"] = "MGH107"
annotation[,c(1)][annotation[,c(1)] == "mgh103"] = "MGH103"

clinical = read.table("/home/wsa/data/young_older_clinical/group7/group7.txt",header = T,sep = "\t",fileEncoding = 'utf-8')
col = colnames(clinical)
sample = clinical$sample
for (j in sample){
  anno_rownumber <- which(annotation$sample == j, arr.ind = TRUE)
  cli_rownumber <- which(clinical$sample == j, arr.ind = TRUE)
  annotation[anno_rownumber,paste(col)] <- clinical[cli_rownumber,]
}


#添加年龄分级信息
anno_rownumber <- which(annotation$age <= 50 &annotation$age >=18, arr.ind = TRUE)
annotation[anno_rownumber,'level'] <- 'adult'
anno_rownumber <- which(annotation$age <= 100 & annotation$age >=65, arr.ind = TRUE)
annotation[anno_rownumber,'level'] <- 'aged'
View(annotation)

#将矩阵和metadata变成seuratobject
scRNAlist_tmp[[7]] <- CreateSeuratObject(counts,min.cells=3, min.features = 200)
scRNAlist_tmp[[7]] <- AddMetaData(
  object = scRNAlist_tmp[[7]],
  metadata = annotation,
  col.name = names(annotation))

#scRNAlist_tmp[[1]]是所有数据，下面选出我想要的数据
tmp <- subset(scRNAlist_tmp[[7]]@meta.data, tumor %in% c("GBM") & level %in% c("adult","aged") & new %in% c("Primary","primary"))
scRNAlist[[7]] <- subset(scRNAlist_tmp[[7]], cells=row.names(tmp))
scRNAlist[[7]] <- RenameCells(scRNAlist[[7]], add.cell.id = 7)
ncol(scRNAlist[[7]])
rm(tmp)







if (FALSE){
  ################################################################################
  #读取GSE151506，多个txt文件
  ##还要修改！！！
  ################################################################################
  #读取矩阵
  t <- list()
  setwd('/home/wsa/data/young_older/GSE151506')
  assays <- dir('/home/wsa/data/young_older/GSE151506')
  dir <- paste0("/home/wsa/data/young_older/GSE151506/", assays)
  
  counts = read.table(dir[[1]], header=TRUE, row.names=1,sep = '\t',check.names = FALSE)
  for (i in 2:length(dir)){
    new.data = read.table(dir[[i]], header=TRUE, row.names=1,sep = '\t',check.names = FALSE)
    counts = cbind(counts,new.data)
  }
  
  
  #读取临床信息
  annotation = read.table("/home/wsa/data/young_older_clinical/GSE151506/anno.txt",header = T,sep='\t',row.names = 1)
  clinical = read.table("/home/wsa/data/young_older_clinical/GSE151506/clinical.txt",header = T,sep = "\t",fileEncoding = 'utf-8')
  col = colnames(clinical)
  sample = clinical$sample
  for (j in sample){
    anno_rownumber <- which(annotation$sample == j, arr.ind = TRUE)
    cli_rownumber <- which(clinical$sample == j, arr.ind = TRUE)
    annotation[anno_rownumber,paste(col)] <- clinical[cli_rownumber,]
  }
  
  #将矩阵和metadata变成seuratobject
  scRNAlist_tmp[[6]] <- CreateSeuratObject(counts,min.cells=3, min.features = 200)
  scRNAlist_tmp[[6]] <- AddMetaData(
    object = scRNAlist_tmp[[6]],
    metadata = annotation,
    col.name = names(annotation))
  
  #scRNAlist_tmp[[1]]是所有数据，下面选出我想要的数据
  tmp1 <- subset(scRNAlist_tmp[[6]]@meta.data, age<=50)
  tmp2 <- subset(scRNAlist_tmp[[6]]@meta.data, age<=70)
  tmp2 <- subset(tmp2, age>=60)
  scRNAlist[[6]] <- subset(scRNAlist_tmp[[6]], cells=c(row.names(tmp1),row.names(tmp2)))
  scRNAlist[[6]] <- RenameCells(scRNAlist[[6]], add.cell.id = 6)
  
}


################################################################################
#添加批次meta标签
################################################################################
for(i in 1:length(scRNAlist)){
  scRNAlist[[i]] <- AddMetaData(
    object = scRNAlist[[i]],
    metadata = c(rep(i,scRNAlist[[i]]@assays$RNA@counts@Dim[2])),
    col.name = "group"
  )
}

scRNAlist[[3]] <- AddMetaData(
  object = scRNAlist[[3]],
  metadata = c(rep(3,scRNAlist[[3]]@assays$RNA@counts@Dim[2])),
  col.name = "group")

################################################################################
#合并数据集
################################################################################
#合并多个数据集
#scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
#View(scRNA)

#################################################################################
#分课题组，保存数据为h5Seurat格式
#################################################################################
setwd('/home/wsa/data/scRNAlist_1850_6590_h5seurat')
for(i in 1:length(scRNAlist)){
  filname <- paste0("scRNAlist",i,".h5Seurat")
  SaveH5Seurat(scRNAlist[[i]],filename = filname)
}

filname <- paste0("scRNAlist",3,".h5Seurat")
SaveH5Seurat(scRNAlist[[3]],filename = filname)


if (FALSE){
################################################################################
#数据质控等数据预处理 https://rotterlt.top/4924.html
#https://www.jianshu.com/p/159df5c4ff21
#https://blog.csdn.net/leianuo123/article/details/111499875
#https://www.jianshu.com/p/d4df73ff6a2e
################################################################################
table(scRNA$orig.ident)
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)


################################################################################
#test

cli_data = read.table("/home/wsa/data/2_test.txt",header = T,sep = "\t",fileEncoding = 'utf-8')

scRNAlist[[1]] <- AddMetaData(
      object = scRNAlist[[1]],
      metadata = c(rep(t[1,],scRNAlist[[1]]@assays$RNA@counts@Dim[2])),
      col.name = colnames(t)
    )
z = rep(t[1,],scRNAlist[[1]]@assays$RNA@counts@Dim[2])
Time = rep(time[1],scRNAlist[[1]]@assays$RNA@counts@Dim[2])
View(scRNAlist[[1]]@meta.data)

df1<-data.frame(col1=c(1,2),col2=c(2,3))
df2<-data.frame(col1=c(1,4),col2=c(2,100))
df <- merge(df1, df2, all=TRUE)
View(df)
4555+11683+5950


m <- read.csv('/home/wsa/data/w_anno/GSM4972210_annot.Human.GBM.R1_2_3_4_4nc.csv', header = T, row.names = 1)
# csv矩阵转换成数据框
m = data.frame(m)
n = read.table("/home/wsa/data/3_test.txt",header = T,sep = "\t",fileEncoding = 'utf-8')
z = colnames(n)
h = n$Name
for (j in h){
  v <- which(m$sample==j, arr.ind = TRUE)
  vv <-which(n$Name==j, arr.ind = TRUE)
  m[v,paste(z)] <- n[vv,]
}
m <- subset(m, select = -c(sample))

##########################
# 按meta.data提取数据子集
##########################
zzz <- subset(scRNA@meta.data, Name=="ND5")
scRNAsub <- subset(scRNA, cells=row.names(zzz))

}

