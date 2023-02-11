library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggpubr)
#install.packages("ggrepel")
#install.packages("pheatmap")
library(ggrepel)
library(pheatmap)
library(tidyr)

scRNA <- readRDS("scRNA_young_older.Rds")

scRNAtmp <- subset(scRNA@meta.data, group == 4)
scRNApeople <- subset(scRNA, cells=row.names(scRNAtmp))

# ############################################################################################################
# ## 使用harmony整合数据
# ############################################################################################################

scRNApeople <- NormalizeData(scRNApeople, normalization.method = "LogNormalize", scale.factor = 10000)  # Checkpoint
scRNApeople <- FindVariableFeatures(scRNApeople, selection.method = "vst", nfeatures = 4000)
all.genes <- rownames(scRNApeople)
scRNApeople <- ScaleData(scRNApeople, features = all.genes)
scRNApeople <- RunPCA(scRNApeople, features = VariableFeatures(object = scRNApeople))
# 未经校正的PC中的数据集之间存在明显差异：
DimPlot(scRNApeople, reduction = "pca",group.by = "group",label = T)
# harmony
scRNApeople <- RunHarmony(scRNApeople,group.by.vars = "group",max.iter.harmony = 20 )
DimPlot(scRNApeople, reduction = "harmony",group.by = "group",label = T,pt.size = 0.1)


# 降维聚类
ElbowPlot(scRNApeople, ndims= 50)
pc.num = 1:35
scRNApeople <- RunUMAP(scRNApeople, reduction = "harmony", dims = pc.num)
scRNApeople <- FindNeighbors(scRNApeople, dims = pc.num)  # louvain cluster, graph based//基于PCA的结果
scRNApeople <- FindClusters(scRNApeople, resolution = 0.5)
identity(scRNApeople)
DimPlot(scRNApeople,reduction = "umap",group.by = "seurat_clusters", pt.size = 0.1)
DimPlot(scRNApeople,reduction = "umap",group.by = "seurat_clusters",split.by = "level",pt.size = 0.1)
DimPlot(scRNApeople,reduction = "umap",group.by = "group",split.by =  "group",pt.size = 0.1)

#GBMs为IDH-具有常见抑癌基因和癌基因突变的野生型
FeaturePlot(scRNApeople, features = c('TP53', 'PTEN', 'TERT', 'CDKN2A', 'CDK4','NF1','CD68'))
#marcophage,M1,M2型
FeaturePlot(scRNApeople, features = c('CD68'),split.by = "level")
FeaturePlot(scRNApeople, features = c('CD86', 'CD80'),split.by = "level")
FeaturePlot(scRNApeople, features = c('CD163', 'CD206'),split.by = "level")
#microglia
FeaturePlot(scRNApeople, features = c('CD74', 'LPAR6', 'DOCK8','AIF1' ),split.by = "level")
# T cell
FeaturePlot(scRNApeople, features = c('CD3','PTPRC','CD3E'),split.by = "level",label =T)
# MDSC
FeaturePlot(scRNApeople, features = c('CD15','CD14','HLA-DR'),split.by = "level")
# neutrophils
FeaturePlot(scRNApeople, features = c('CXCL3','CSF3R','S100A8','S100A9'),split.by = "level")
#CD4-CD8+T
FeaturePlot(scRNApeople, features = c('NKG7','CD8A ','Cd3d '))
# Bcell
FeaturePlot(scRNApeople, features = c("Ptprc","Cd3e","Cd4", "Cd8a", "Cd79a", "Cd79b",  "Nkg7", "Itgam", "Csf1r","Itgax","Flt3","Csf3r"), reduction = "umap",cols = c("gray", "red"))
# 免疫细胞
FeaturePlot(scRNApeople, features = c('PTPRC'),split.by = 'group')
# tumor
FeaturePlot(scRNApeople, features = c('PTPRZ1','OLIG2 ','AQP4','HYDIN','MKI67'),split.by = 'group')


# ############################################################################################################
# ## QC
# ############################################################################################################
#计算线粒体基因比例
#scRNApeople[["percent.mt"]] <- PercentageFeatureSet(scRNApeople, pattern = "^MT-")
#VlnPlot(scRNApeople, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#plot1 <- FeatureScatter(scRNApeople, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(scRNApeople, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2

#scRNApeople <- subset(scRNApeople, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# ############################################################################################################
# ## Normalization
# ############################################################################################################
scRNApeople <- NormalizeData(scRNApeople, normalization.method = "LogNormalize", scale.factor = 10000)  # Checkpoint

# ############################################################################################################
# ## Feature selection
# ############################################################################################################
scRNApeople <- FindVariableFeatures(scRNApeople, selection.method = "vst", nfeatures = 4000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scRNApeople), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(scRNApeople)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

all.genes <- rownames(scRNApeople)
scRNApeople <- ScaleData(scRNApeople, features = all.genes)
# ############################################################################################################
# ## Reduction
# ############################################################################################################
scRNApeople <- RunPCA(scRNApeople, features = VariableFeatures(object = scRNApeople))
DimPlot(scRNApeople, reduction = "pca")
# FeaturePlot(scRNApeople, reduction = 'pca', features = c('CD79A', 'CD14', 'FCGR3A', 'CD4', 'CST3', 'PPBP'))
# FeaturePlot(scRNApeople, reduction = 'pca', features = c('nCount_RNA', 'nFeature_RNA'))
# loadings_sorted = dplyr::arrange(as.data.frame(loadings), desc(PC_1))
DimHeatmap(scRNApeople, dims = 1:2, cells = 200, balanced = TRUE)
ElbowPlot(scRNApeople, ndims = 50)
scRNApeople <- RunUMAP(scRNApeople, dims = 1:25) # umap tsne
FeaturePlot(scRNApeople, features = c('FCGR3A', 'CD14'), reduction = 'umap')
# ############################################################################################################
# ## Cluster
# ############################################################################################################
scRNApeople <- FindNeighbors(scRNApeople, dims = 1:25)  # louvain cluster, graph based//基于PCA的结果
scRNApeople <- FindClusters(scRNApeople, resolution = .5)
DimPlot(scRNApeople, reduction = "umap",group.by = "group",label = T)
DimPlot(scRNApeople, reduction = "umap", group.by = 'seurat_clusters', label=T)

# ############################################################################################################
# ## 注释
# ############################################################################################################
FeaturePlot(scRNApeople, features = c("MS4A1", "TYROBP", "CD14",'FCGR3A', "FCER1A",
                                "CCR7", "IL7R", "PPBP", "CD8A"))
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(scRNApeople)
scRNApeople <- RenameIdents(scRNApeople, new.cluster.ids)
DimPlot(scRNApeople, reduction = "umap", label=T)


################################################################################
## 加载seurat，将umap坐标加入到seurat对象的reducObject
################################################################################

seurat <- SeuratDisk::LoadH5Seurat("/home/wsa/data/all_1850_6590.h5seurat")

alt_cell_type <- read.table('/home/wsa/data/mtx_file_1850_6590/alt_cluster_predict.txt',sep = ',', header = T, row.names = 1,col.names = "V1")
ref_cell_type <- read.table('/home/wsa/data/mtx_file_1850_6590/ref_dataset_cluster.txt',sep = '\n')
cell_type = rbind(alt_cell_type,ref_cell_type)
rownames(cell_type) <- colnames(seurat)
nrow(cell_type)
nrow(seurat@meta.data)
seurat <- AddMetaData(
  object = seurat,
  metadata = cell_type,
  col.name = "cell_type")

umap1 <- read.table('/home/wsa/data/bench_glioma_1850_6590_umap.txt',sep=',', header = TRUE,row.names = 1)
rownames(umap1) = colnames(seurat)
colnames(umap1) <- paste0("UMAP1_", 1:2)
seurat[["umap1"]] <- CreateDimReducObject(embeddings = as.matrix(umap1), key = "UMAP1_",assay = DefaultAssay(seurat))
DimPlot(seurat,reduction="umap1",group.by = "cell_type",pt.size = 0.1)
DimPlot(seurat,reduction="umap1",group.by = "cell_type",split.by = "group",pt.size = 0.1)
theme_feature <- theme(axis.line = element_blank(), 
                       axis.text = element_blank(), 
                       axis.ticks = element_blank(),
                       axis.title = element_blank(),
                       plot.title = element_text(size = 20,hjust = 0.5, face="bold",colour = "#000000")) 
p1 <- FeaturePlot(seurat, features = c('EGFR'),reduction = 'umap1', max.cutoff = 3, cols = c("#f2e9e4","#F82547"),slot = "data",label.size = 6,pt.size = 0.1)+ 
  theme_feature+
  ggtitle("EGFR")
p2 <- FeaturePlot(seurat, features = c('PTPRC'),reduction = 'umap1', max.cutoff = 3, cols = c("#f2e9e4","#F82547"),slot = "data",label.size = 6,pt.size = 0.1)+
  theme_feature+
  ggtitle("PTPRC(CD45)")
p <- p1 + p2 + #加号必须在第一行表示代码还没结束 
  plot_annotation(tag_levels = "A")+
  plot_layout(ncol = 2,#图形设置为两列，默认按行填充，
              widths = c(1, 1))#两列之间相对宽度比为2：1

FeaturePlot(seurat, features = c('IL2RA'),reduction = 'umap1', max.cutoff = 5, cols = c("#f2e9e4","#F82547"),slot = "data",label.size = 6,pt.size = 0.1)+ 
  theme_feature
FeaturePlot(seurat, features = c('SCIN'),reduction = 'umap1', max.cutoff = 5, cols = c("#F9F871","#FFB349", "#F76D4D","#CF2561","#8A0072","#020f75"),slot = "data",label.size = 6,pt.size = 0.1)+ 
  theme_feature


################################################################################
## 细胞再聚类
################################################################################
Cells.sub <- subset(seurat@meta.data, cell_type=="TAM 2")
scRNAsub <- subset(seurat, cells=row.names(Cells.sub))
##PCA降维
scRNAsub <- FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 2000)
scale.genes <-  rownames(scRNAsub)
scRNAsub <- ScaleData(scRNAsub, features = scale.genes)
scRNAsub <- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub))
ElbowPlot(scRNAsub, ndims=50,reduction = "pca")
pc.num=1:25
scRNAsub <- RunHarmony(scRNAsub,group.by.vars = "group",max.iter.harmony = 40, lambda = 0.2, theta = 0.5 ,dim.use = pc.num, sigma = 1)

##非线性降维
#UMAP
scRNAsub <- RunUMAP(scRNAsub, dims = pc.num, reduction = "harmony")
scRNAsub <- FindNeighbors(scRNAsub, dims = pc.num, reduction = "harmony") 
scRNAsub <- FindClusters(scRNAsub, resolution = 0.2)
table(scRNAsub@meta.data$seurat_clusters)
#metadata <- scRNAsub@meta.data
#cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
#write.csv(cell_cluster,'subcluster/cell_cluster.csv',row.names = F)

embed_umap <- Embeddings(scRNAsub, 'umap')
#write.csv(embed_umap,'subcluster/embed_umap.csv') 
DimPlot(scRNAsub, reduction = "umap",pt.size = 0.1) 
DimPlot(scRNAsub, reduction = "umap",group.by = "seurat_clusters",split.by = "group")

scRNA_singleR <- GetAssayData(scRNAsub, slot="data")
hesc <- SingleR(test = scRNA_singleR, ref = list(ref_database_immune,ref_monaco_immune), labels = list(ref_database_immune$label.fine,ref_monaco_immune$label.fine)) 
hesc
#seurat 和 singleR的table表
table(hesc$labels)
scRNAsub@meta.data$l <- hesc$labels
DimPlot(scRNAsub, group.by = c("cell_type", "l"),reduction = "umap")

ggsave("subcluster/UMAP.pdf", plot = plot2, width = 8, height = 7)
ggsave("subcluster/UMAP.png", plot = plot2, width = 8, height = 7)



# ##############################################################################
# ## 寻找差异基因，做火山图
# ##############################################################################
#修改ident名字为细胞类型
levels(seurat)  #查看当前Idents
Idents(seurat)=seurat@meta.data$cell_type  #这句代码切换
levels(seurat)  #查看是否更改
seurat@meta.data$orig.ident <- NULL
seurat@meta.data$orig.ident <- seurat@meta.data$cell_type


###################
## 不同类型细胞比较
markers <- FindMarkers(seurat, ident.1 = "TAM 1", ident.2 = "Tumor", min.pct = 0.25,logfc.threshold = 0.1,only.pos = FALSE)
#setwd('/home/wsa/data/differential_gene/')
#write.table(markers,"MDM_DG.csv",row.names=TRUE,col.names=TRUE,sep=",")
markers$Group = "NoSignifi"
markers$Group[which((markers$p_val_adj < 0.001) & (markers$avg_log2FC > 0.4))] = "Up"
markers$Group[which((markers$p_val_adj < 0.001) & (markers$avg_log2FC < -0.4))] = "Down"
table(markers$Group)
data <- markers
data$gene <- row.names(data)
data$logP <- -log10(data$p_val_adj)
ggplot(data = data, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Group, label = gene)) + 
  geom_point(alpha=0.8, size=4)+
  scale_color_manual(values=c("#4393C3","#00000033","#FC4E2A"))+
  ylim(c(-1, 310))+
  labs(x="log2(fold change)",y="-log10 (p-value)",title="microglia vs MDM")+
  geom_text_repel(
    data = subset(data, data$logP > 40 & abs(data$avg_log2FC) >= 0.8),
    aes(label = gene),
    size = 3,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), 
    color = "black",
    segment.color = "black", show.legend = FALSE )+
  vol_theme

mark <- FindMarkers(seurat, ident.1 = "TAM 1", ident.2 = "Tumor",min.pct = 0.25,logfc.threshold = 0.4,only.pos = FALSE)
View(mark)
setwd('/home/wsa/data/differential_gene/')
write.table(mark,"MDM_tumor.csv",row.names=TRUE,col.names=TRUE,sep=",")

mark <- FindMarkers(seurat, ident.1 = "TAM 2", ident.2 = "Tumor",min.pct = 0.25,logfc.threshold = 0.4,only.pos = FALSE)
mark
setwd('/home/wsa/data/differential_gene/')
write.table(mark,"microglia_tumor.csv",row.names=TRUE,col.names=TRUE,sep=",")

mark <- FindMarkers(seurat, ident.1 = "T cells", ident.2 = "Tumor",min.pct = 0.25,logfc.threshold = 0.4,only.pos = FALSE)
mark
setwd('/home/wsa/data/differential_gene/')
write.table(mark,"Tcell_tumor.csv",row.names=TRUE,col.names=TRUE,sep=",")

mark <- FindMarkers(seurat, ident.1 = "NK cells", ident.2 = "Tumor",min.pct = 0.25,logfc.threshold = 0.4,only.pos = FALSE)
mark
setwd('/home/wsa/data/differential_gene/')
write.table(mark,"NKcell_tumor.csv",row.names=TRUE,col.names=TRUE,sep=",")



#######################
## 同一细胞不同年龄比较
vol_theme <-   theme(panel.background = element_blank(),  #去除灰色背景
                     panel.grid = element_blank(), #去除网格线
                     axis.line = element_line(colour = "#000000",size = 0.3),
                     axis.text = element_text(colour = "#000000" ,size = 18),
                     #axis.text.x = element_text(angle = 60,vjust = 1,hjust = 1), 
                     axis.ticks = element_line(colour = "#000000" ,size = 0.3),
                     axis.ticks.length = unit(1,'mm'),
                     plot.margin = unit(c(0.5,0.4,0.4,0.3),"cm"),
                     axis.title.y = element_text(size = 20),
                     axis.title.x = element_blank(),
                     plot.title = element_text(size = 24,hjust = 0.5),
                     legend.text = element_text(size=16), #设置legend标签的大小
                     legend.key.size=unit(0.8,'cm')  )  # 设置legend标签之间的大小


Cells.sub1 <- subset(seurat@meta.data, cell_type %in% "TAM 1")
scRNAsub1 <- subset(seurat, cells=row.names(Cells.sub1))
levels(scRNAsub1)  #查看当前Idents
Idents(scRNAsub1)=scRNAsub1@meta.data$level  #这句代码切换
levels(scRNAsub1)  #查看是否更改
scRNAsub1@meta.data$orig.ident <- NULL
scRNAsub1@meta.data$orig.ident <- scRNAsub1@meta.data$level
TAM1.markers <- FindMarkers(scRNAsub1, ident.1 = "adult", ident.2 = "aged",min.pct = 0.25,logfc.threshold = 0.1,only.pos = FALSE)
View(TAM1.markers)
#setwd('/home/wsa/data/differential_gene/')
#write.table(TAM1.markers,"MDM_young_older_DG.csv",row.names=TRUE,col.names=TRUE,sep=",")

TAM1.markers$Group = "NoSignifi"
TAM1.markers$Group[which((TAM1.markers$p_val_adj < 0.001) & (TAM1.markers$avg_log2FC > 0.4))] = "Up"
TAM1.markers$Group[which((TAM1.markers$p_val_adj < 0.001) & (TAM1.markers$avg_log2FC < -0.4))] = "Down"
table(TAM1.markers$Group)
data <- TAM1.markers
data$gene <- row.names(data)
data$logP <- -log10(data$p_val_adj)

ggplot(data = data, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Group, label = gene)) + 
  geom_point(alpha=0.8, size=4)+
  scale_color_manual(values=c("#4393C3","#00000033","#FC4E2A"))+
  ylim(c(-1, 310))+
  labs(x="log2(fold change)",y="-log10 (p-value)",title="MDM aged vs adult")+
  geom_text_repel(
    data = subset(data, data$logP > 10 & abs(data$avg_log2FC) >= 0.7),
    aes(label = gene),
    size = 3,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), 
    color = "black",
    max.overlaps = 20,
    segment.color = "black", show.legend = FALSE )+
  vol_theme


Cells.sub2 <- subset(seurat@meta.data, cell_type %in% "TAM 2")
scRNAsub2 <- subset(seurat, cells=row.names(Cells.sub2))
table(scRNAsub2@meta.data$level)
levels(scRNAsub2)  #查看当前Idents
Idents(scRNAsub2)=scRNAsub2@meta.data$level  #这句代码切换
levels(scRNAsub2)  #查看是否更改
scRNAsub2@meta.data$orig.ident <- NULL
scRNAsub2@meta.data$orig.ident <- scRNAsub2@meta.data$level
TAM2.markers <- FindMarkers(scRNAsub2, ident.1 = "adult", ident.2 = "aged",min.pct = 0.25,logfc.threshold = 0.1,only.pos = FALSE)
View(TAM2.markers)
#setwd('/home/wsa/data/differential_gene/')
#write.table(TAM2_young.markers,"microglia_young_DG.csv",row.names=TRUE,col.names=TRUE,sep=",")

TAM2.markers$Group = "NoSignifi"
TAM2.markers$Group[which((TAM2.markers$p_val_adj < 0.001) & (TAM2.markers$avg_log2FC > 0.4))] = "Up"
TAM2.markers$Group[which((TAM2.markers$p_val_adj < 0.001) & (TAM2.markers$avg_log2FC < -0.4))] = "Down"
table(TAM2.markers$Group)
data <- TAM2.markers
data$gene <- row.names(data)
data$logP <- -log10(data$p_val_adj)
ggplot(data = data, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Group, label = gene)) + 
  geom_point(alpha=0.8, size=4)+
  scale_color_manual(values=c("#4393C3","#00000033","#FC4E2A"))+
  ylim(c(-1, 310))+
  labs(x="log2(fold change)",y="-log10 (p-value)",title="microglia aged vs adult")+
  geom_text_repel(
    data = subset(data, data$logP > 40 & abs(data$avg_log2FC) >= 0.6),
    aes(label = gene),
    size = 3,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), 
    color = "black",
    max.overlaps = 20,
    segment.color = "black", show.legend = FALSE )+
  vol_theme



Cells.sub1 <- subset(seurat@meta.data, cell_type %in% "Tumor")
scRNAsub1 <- subset(seurat, cells=row.names(Cells.sub1))
levels(scRNAsub1)  #查看当前Idents
Idents(scRNAsub1)=scRNAsub1@meta.data$level  #这句代码切换
levels(scRNAsub1)  #查看是否更改
scRNAsub1@meta.data$orig.ident <- NULL
scRNAsub1@meta.data$orig.ident <- scRNAsub1@meta.data$level
markers <- FindMarkers(scRNAsub1, ident.1 = "adult", ident.2 = "aged",min.pct = 0.25,logfc.threshold = 0.4,only.pos = FALSE)
View(markers)
setwd('/home/wsa/data/differential_gene/')
write.table(markers,"Tumor_adult_aged.csv",row.names=TRUE,col.names=TRUE,sep=",")

Cells.sub1 <- subset(seurat@meta.data, cell_type %in% "T cells")
scRNAsub1 <- subset(seurat, cells=row.names(Cells.sub1))
levels(scRNAsub1)  #查看当前Idents
Idents(scRNAsub1)=scRNAsub1@meta.data$level  #这句代码切换
levels(scRNAsub1)  #查看是否更改
scRNAsub1@meta.data$orig.ident <- NULL
scRNAsub1@meta.data$orig.ident <- scRNAsub1@meta.data$level
markers <- FindMarkers(scRNAsub1, ident.1 = "adult", ident.2 = "aged",min.pct = 0.25,logfc.threshold = 0.4,only.pos = FALSE)
markers
setwd('/home/wsa/data/differential_gene/')
write.table(markers,"Tcell_adult_aged.csv",row.names=TRUE,col.names=TRUE,sep=",")

Cells.sub1 <- subset(seurat@meta.data, cell_type %in% "NK cells")
scRNAsub1 <- subset(seurat, cells=row.names(Cells.sub1))
levels(scRNAsub1)  #查看当前Idents
Idents(scRNAsub1)=scRNAsub1@meta.data$level  #这句代码切换
levels(scRNAsub1)  #查看是否更改
scRNAsub1@meta.data$orig.ident <- NULL
scRNAsub1@meta.data$orig.ident <- scRNAsub1@meta.data$level
markers <- FindMarkers(scRNAsub1, ident.1 = "adult", ident.2 = "aged",min.pct = 0.25,logfc.threshold = 0.4,only.pos = FALSE)
markers
setwd('/home/wsa/data/differential_gene/')
write.table(markers,"NKcell_adult_aged.csv",row.names=TRUE,col.names=TRUE,sep=",")

Cells.sub1 <- subset(seurat@meta.data, cell_type %in% "TAM 1")
scRNAsub1 <- subset(seurat, cells=row.names(Cells.sub1))
levels(scRNAsub1)  #查看当前Idents
Idents(scRNAsub1)=scRNAsub1@meta.data$level  #这句代码切换
levels(scRNAsub1)  #查看是否更改
scRNAsub1@meta.data$orig.ident <- NULL
scRNAsub1@meta.data$orig.ident <- scRNAsub1@meta.data$level
markers <- FindMarkers(scRNAsub1, ident.1 = "adult", ident.2 = "aged",min.pct = 0.25,logfc.threshold = 0.4,only.pos = FALSE)
markers
setwd('/home/wsa/data/differential_gene/')
write.table(markers,"MDM_adult_aged.csv",row.names=TRUE,col.names=TRUE,sep=",")

Cells.sub1 <- subset(seurat@meta.data, cell_type %in% "TAM 2")
scRNAsub1 <- subset(seurat, cells=row.names(Cells.sub1))
levels(scRNAsub1)  #查看当前Idents
Idents(scRNAsub1)=scRNAsub1@meta.data$level  #这句代码切换
levels(scRNAsub1)  #查看是否更改
scRNAsub1@meta.data$orig.ident <- NULL
scRNAsub1@meta.data$orig.ident <- scRNAsub1@meta.data$level
markers <- FindMarkers(scRNAsub1, ident.1 = "adult", ident.2 = "aged",min.pct = 0.25,logfc.threshold = 0.4,only.pos = FALSE)
markers
setwd('/home/wsa/data/differential_gene/')
write.table(markers,"microglia_adult_aged.csv",row.names=TRUE,col.names=TRUE,sep=",")







FeaturePlot(seurat, features = c('HMOX1','PLIN2'),reduction = 'umap1', max.cutoff = 7) #缺氧
FeaturePlot(seurat, features = c('CXCL2'),reduction = 'umap1', max.cutoff = 3) #Cytokines and inflammatory response
FeaturePlot(seurat, features = c('P2RY12','TMEM119'),reduction = 'umap1', max.cutoff = 3) #小胶质细胞
FeaturePlot(seurat, features = c("CD14","FCGR3A","CDKN1C"),reduction = 'umap1', max.cutoff = 3,split.by = 'level') 
FeaturePlot(seurat, features = c("TMIGD3","APOC2","SCIN"),reduction = 'umap1', max.cutoff = 3,split.by = 'level') 




# ##############################################################################
# ## 做热图
# ##############################################################################
seurat <- ScaleData(object = seurat, features = rownames(seurat))
set.seed(42)
subobj <- subset(seurat, downsample = 1000)
subobj$cell_type <- factor(x=subobj$cell_type,
                          levels = c("Tumor","Monocytes","TAM 1","TAM 2","prol. TAM","NormalBrain","DC","B cells","T cells","NK cells"))
merged.markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4)
top10 <- merged.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

options(repr.plot.width = 13, repr.plot.height=13)
DoHeatmap(subobj, features = top10$gene, disp.min=-2.5, disp.max=2.5)+
  scale_fill_gradientn(colors = c("#00A8CE","#18B4CA","#EEE8A9","red"))

##tumor normal monocyte macrophage DC T NK B
genelist <- c("EGFR","PTN","SOX4","SEC61G","CLU","CDK4","PTPRZ1","FABP7","CHI3L1","C1orf61",    #tumor
              "S100A9","S100A8","LYZ","S100A12","S100A4","AREG","FCN1","EREG","TIMP1","CFP","VCAN","APOBEC3A","SPN","CCL2","CCL7",   #monocyte
              "C1QA","APOC1","APOE","MAFB","CD14","FCGR1A","CSF1R","CD68","CD14","MAFB",   #macrophage
              "PLP1"  ,   "PTGDS"   , "PPP1R14A", "TF"   ,    "MBP" ,     "MAG"  ,    "CNP"   ,   "QDPR"   ,  "APOD",   #normalBrain
              "CLEC9A","TACSTD2","CLNK","IDO1","FCER1A","CD1C","CD1E","CLEC10A","CCR7","CCL22","IL3RA","LILRA4","IRF7",  #DC
              "MS4A1","BANK1","CD79A","IGHM",     "JCHAIN",   "IGHG1",  "IGHA1"  ,  "IGLC3" ,   "IGKC"  ,   "IGLC2","CD79B","CD19",  #B cell
              "CD3G","TRAC","TRBC2",   #T cell
              "KLRC1","TRDC","GNLY", "GZMA","PRF1","NKG7"     #NK cell
              )

DoHeatmap(subobj, features = genelist, disp.min=-2.5, disp.max=2.5,lines.width = 70,draw.lines = TRUE,group.by = "cell_type")+
  scale_fill_gradientn(colors = c("#00A8CE","#18B4CA","#EEE8A9","red"))+ NoLegend()




mat = GetAssayData(subobj, slot="counts")
p <- pheatmap(mat)



# ##############################################################################
# ## 做细胞因子气泡图
# ##############################################################################

seurat@meta.data$cell_type[which(seurat@meta.data$cell_type =='TAM 1')] <- 'MDM'
seurat@meta.data$cell_type[which(seurat@meta.data$cell_type =='TAM 2')] <- 'Microglia'
seurat@meta.data$cell_type[which(seurat@meta.data$cell_type =='DC')] <- 'Dendritic Cells'
seurat@meta.data$cell_type[which(seurat@meta.data$cell_type =='prol. TAM')] <- 'TAM 3'
table(seurat@meta.data$cell_type)
cytokine_marker <- read.table('/home/wsa/data/differential_gene/tumor_cytokine2.txt',header = F,sep = "\n")
cytokine <- tidyr::unite(seurat@meta.data, "cytokine", cell_type,level, remove = FALSE)
cytokine <- cytokine[,"cytokine"]
cytokine <- as.data.frame(cytokine,row.names = rownames(seurat@meta.data))
seurat <- AddMetaData(
  object = seurat,
  metadata = cytokine,
  col.name = "cytokine")
tmp <- subset(seurat@meta.data, !(cytokine %in% c("NormalBrain_NA","Tumor_NA")))
dot <- subset(seurat, cells=row.names(tmp))
DotPlot(dot, features = cytokine_marker$V1,group.by ="cytokine")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_gradientn(colours = c("#00A8CE","#18B4CA","red"))+
  theme(panel.grid.major.x = element_line(colour = "grey",linetype = 1),
        axis.text.x = element_text(size = 8, angle = 60,vjust = 1,hjust = 1),
        axis.text.y = element_text(size = 10,vjust = 1,hjust = 1))


# ##############################################################################
# ## singleR给TAM打细分标签
# ##############################################################################
#library(celldex)
#HumanPrimaryCellAtlasData()
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
#install.packages("remotes")
#remotes::install_github("LTLA/celldex")
#library(SingleR)
#library(celldex)
#library(Seurat)
#library(pheatmap)

#直接load下载好的数据库
setwd('/home/wsa/data/celldex_ref')
load("ref_Monaco_114s.RData")
ref_monaco_immune <- ref_Monaco
ref_database_immune <- readRDS("ref.rds")

Cells_singleR <- subset(seurat@meta.data, cell_type=="TAM 1")
scRNA <- subset(seurat, cells=row.names(Cells_singleR))
scRNA_singleR <- GetAssayData(scRNA, slot="data")
hesc <- SingleR(test = scRNA_singleR, ref = list(ref_database_immune,ref_monaco_immune), labels = list(ref_database_immune$label.fine,ref_monaco_immune$label.fine)) 
hesc
#seurat 和 singleR的table表
table(hesc$labels)

scRNA@meta.data$l <- hesc$labels
DimPlot(scRNA, group.by = c("cell_type", "l"),reduction = "umap1",pt.size = 1)









