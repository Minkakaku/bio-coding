
library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(scales)
library(RColorBrewer)
library(ggsci)

###############################################
## 读取训练好的position
###############################################
mydata <- read.table('/home/wsa/data/bench_glioma_1850_6590_umap.txt',sep=',', header = TRUE,row.names = 1)

alt_cell_type <- read.table('/home/wsa/data/mtx_file_1850_6590/alt_cluster_predict.txt',sep = ',', header = T, row.names = 1,col.names = "V1")
ref_cell_type <- read.table('/home/wsa/data/mtx_file_1850_6590/ref_dataset_cluster.txt',sep = '\n')
cell_type = rbind(alt_cell_type,ref_cell_type)
mydata$cell_type <- cell_type$V1
mydata$cell_type[which(mydata$cell_type =='TAM 1')] <- 'MDM'
mydata$cell_type[which(mydata$cell_type =='TAM 2')] <- 'Microglia'
mydata$cell_type[which(mydata$cell_type =='DC')] <- 'Dendritic Cells'
mydata$cell_type[which(mydata$cell_type =='prol. TAM')] <- 'TAM 3'
table(mydata$cell_type)
level <- c(seurat.1@meta.data$level,seurat.2@meta.data$level,seurat.3@meta.data$level,seurat.5@meta.data$level,seurat.7@meta.data$level,seurat.4@meta.data$level,rep("NA",ncol(seurat.6)))
mydata$level <- level
sample <- c(seurat.1@meta.data$sample,seurat.2@meta.data$sample,seurat.3@meta.data$sample,seurat.5@meta.data$sample,seurat.7@meta.data$sample,seurat.4@meta.data$sample,seurat.6@meta.data$sample)
mydata$sample <- sample
age <- seurat@meta.data$age
mydata$age <- age
group <- seurat@meta.data$group
mydata$group <- group
mydata$cell_type2 <- mydata$cell_type
mydata[which(!(mydata$cell_type2 %in% c("NormalBrain","Tumor"))),'cell_type2'] <- 'immune'

###############################################
## 挑选出TAM1的加入sublabel
###############################################
subumap_data <- subset(mydata, cell_type %in% c('prol. TAM'))
setwd('/home/wsa/data/subembed/')
write.table(subumap_data[1:2], file = 'subpos_TAM3.txt',sep = ',', row.names = F, col.names = F)
sublabel <- read.table('/home/wsa/data/sublabel/sublabel_TAM3.txt',sep = ',', header = T, row.names = 1)
subumap_data$sublabel <- sublabel[,1]

##############################################
## 读取训练好的embed,并挑选出TAM1用于python亚聚类
##############################################
embed_data <- read.table('/home/wsa/data/bench_glioma_with_tumor_super_dataset_scp_embed.txt',sep=',', header = TRUE,row.names = 1) 
alt_cell_type <- read.table('/home/wsa/data/mtx_file/alt_cluster_predict.txt',sep = ',', header = T, row.names = 1,col.names = "V1")
ref_cell_type <- read.table('/home/wsa/data/mtx_file/ref_dataset_cluster.txt',sep = '\n')
cell_type = rbind(alt_cell_type,ref_cell_type)
embed_data$cell_type <- cell_type$V1
View(embed_data)
subembed <- subset(embed_data, cell_type %in% c('prol. TAM'))
setwd('/home/wsa/data/subembed/')
write.table(subembed[1:50], file = 'subembed_TAM3.txt',sep = ',', row.names = F, col.names = F)


################################################################################
## 聚类图
################################################################################
## 所有年龄层
cluster_theme <- theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"),
        legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=18), #设置legend标签的大小
        legend.key.size=unit(1,'cm')  ,  # 设置legend标签之间的大小
        plot.title = element_text(size = 24,hjust = 0.5,colour = "#000000"))
  

p1 <- ggplot2::ggplot(mydata,aes(x = X0, y = X1, color = cell_type))+
  geom_point(size = 0.1, alpha = 1) +
  scale_color_manual(values = c("#91D1C2FF","#B09C85FF","#8491B4FF","#DC0000FF","#7E6148FF","#00A087FF","#4DBBD5FF","#F39B7FFF","#3C5488FF","#E64B35FF"),
                     labels=c("TAM 2" = "microglia","TAM 1" = "MDM","prol. TAM" = "TAM 3"))+
  #scale_color_manual(values = c("brown3","chartreuse3","khaki1","plum1","lightskyblue","darkorange","hotpink","seashell3","#E87796","#E64B35FF"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+ #设置legend中 点的大小 
  ggtitle("Our Cell Clusters") +
  cluster_theme +
  #绘制左下角箭头
  geom_segment(aes(x = min(mydata$X0) , y = min(mydata$X1) ,
                   xend = min(mydata$X0) +6, yend = min(mydata$X1) ),
               colour = "black", size=0.3,arrow = arrow(length = unit(0.3,"cm"),type="closed",angle=15))+ 
  geom_segment(aes(x = min(mydata$X0)  , y = min(mydata$X1)  ,
                   xend = min(mydata$X0) , yend = min(mydata$X1) + 6),
               colour = "black", size=0.3,arrow = arrow(length = unit(0.3,"cm"),type="closed",angle=15)) +
  annotate("text", x = min(mydata$X0) +2, y = min(mydata$X1) -1, label = "UMAP_1",
           color="black",size = 5 ) + 
  annotate("text", x = min(mydata$X0) -1, y = min(mydata$X1) + 2, label = "UMAP_2",
           color="black",size = 5 ,angle=90)

p2 <- ggplot2::ggplot(mydata,aes(x = X0, y = X1, color = factor(group)))+
  geom_point(size = 0.1, alpha = 1) +
  scale_color_manual(values = c("#91D1C2FF","#B09C85FF","#8491B4FF","#F39B7FFF","#7E6148FF","#00A087FF","#4DBBD5FF","#3C5488FF","#E64B35FF"),
                     labels=c("TAM 2" = "microglia","TAM 1" = "MDM","prol. TAM" = "TAM 3"))+
  #scale_color_manual(values = c("brown3","chartreuse3","khaki1","plum1","lightskyblue","darkorange","hotpink","seashell3","#E87796","#E64B35FF"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+ #设置legend中 点的大小 
  ggtitle("Our Group Clusters") +
  cluster_theme +
  #绘制左下角箭头
  geom_segment(aes(x = min(mydata$X0) , y = min(mydata$X1) ,
                   xend = min(mydata$X0) +6, yend = min(mydata$X1) ),
               colour = "black", size=0.3,arrow = arrow(length = unit(0.3,"cm"),type="closed",angle=15))+ 
  geom_segment(aes(x = min(mydata$X0)  , y = min(mydata$X1)  ,
                   xend = min(mydata$X0) , yend = min(mydata$X1) + 6),
               colour = "black", size=0.3,arrow = arrow(length = unit(0.3,"cm"),type="closed",angle=15)) +
  annotate("text", x = min(mydata$X0) +2, y = min(mydata$X1) -1, label = "UMAP_1",
           color="black",size = 5 ) + 
  annotate("text", x = min(mydata$X0) -1, y = min(mydata$X1) + 2, label = "UMAP_2",
           color="black",size = 5 ,angle=90)


setwd('/home/wsa/result/')
pdf('Our Cell Clusters.pdf',width = 8,height = 6)
p1
dev.off()
pdf('Our Group Clusters.pdf',width = 8,height = 6)
p2
dev.off()

png("Our Cell Clusters.png",units="in", width=8, height=6,res=300)
p1
dev.off()
png("Our Group Clusters.png",units="in", width=8, height=6,res=300)
p2
dev.off()


## 免疫细胞与肿瘤细胞聚类图--年轻人
group <- seurat@meta.data$group
mydata$group <- group
mydata$cell_type2 <- mydata$cell_type
mydata[which(!(mydata$cell_type2 %in% c("NormalBrain","Tumor"))),'cell_type2'] <- 'immune'
young_data <- subset(mydata,!(group %in% c('4','6')))
young_data<-subset(young_data, level %in% 'adult')
nrow(young_data)

p1 <- ggplot2::ggplot(young_data,aes(x = X0, y = X1, color = cell_type))+
  geom_point(size = 0.1, alpha = 1) +
  scale_color_manual(values = c("#91D1C2FF","#B09C85FF","#8491B4FF","#DC0000FF","#7E6148FF","#00A087FF","#4DBBD5FF","#F39B7FFF","#3C5488FF","#E64B35FF"),
                     labels=c("TAM 2" = "microglia","TAM 1" = "MDM","prol. TAM" = "TAM 3"))+
  #scale_color_manual(values = c("brown3","chartreuse3","khaki1","plum1","lightskyblue","darkorange","hotpink","seashell3","#E87796","#E64B35FF"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+ #设置legend中 点的大小 
  ggtitle("Cells in Adult") +
  cluster_theme +
  #绘制左下角箭头
  geom_segment(aes(x = min(mydata$X0) , y = min(mydata$X1) ,
                   xend = min(mydata$X0) +6, yend = min(mydata$X1) ),
               colour = "black", size=0.3,arrow = arrow(length = unit(0.3,"cm"),type="closed",angle=15))+ 
  geom_segment(aes(x = min(mydata$X0)  , y = min(mydata$X1)  ,
                   xend = min(mydata$X0) , yend = min(mydata$X1) + 6),
               colour = "black", size=0.3,arrow = arrow(length = unit(0.3,"cm"),type="closed",angle=15)) +
  annotate("text", x = min(mydata$X0) +2, y = min(mydata$X1) -1, label = "UMAP_1",
           color="black",size = 5 ) + 
  annotate("text", x = min(mydata$X0) -1, y = min(mydata$X1) + 2, label = "UMAP_2",
           color="black",size = 5 ,angle=90)

## 免疫细胞与肿瘤细胞聚类图--老年人
older_data <- subset(mydata,!(group %in% c('4','6')))
older_data<-subset(older_data, level %in% 'aged')
nrow(older_data)
p2 <- ggplot2::ggplot(older_data,aes(x = X0, y = X1, color = cell_type))+
  geom_point(size = 0.1, alpha = 1) +
  scale_color_manual(values = c("#91D1C2FF","#B09C85FF","#8491B4FF","#DC0000FF","#7E6148FF","#00A087FF","#4DBBD5FF","#F39B7FFF","#3C5488FF","#E64B35FF"),
                     labels=c("TAM 2" = "microglia","TAM 1" = "MDM","prol. TAM" = "TAM 3"))+
  #scale_color_manual(values = c("brown3","chartreuse3","khaki1","plum1","lightskyblue","darkorange","hotpink","seashell3","#E87796","#E64B35FF"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+ #设置legend中 点的大小 
  ggtitle("Cells in Aged") +
  cluster_theme +
  #绘制左下角箭头
  geom_segment(aes(x = min(mydata$X0) , y = min(mydata$X1) ,
                   xend = min(mydata$X0) +6, yend = min(mydata$X1) ),
               colour = "black", size=0.3,arrow = arrow(length = unit(0.3,"cm"),type="closed",angle=15))+ 
  geom_segment(aes(x = min(mydata$X0)  , y = min(mydata$X1)  ,
                   xend = min(mydata$X0) , yend = min(mydata$X1) + 6),
               colour = "black", size=0.3,arrow = arrow(length = unit(0.3,"cm"),type="closed",angle=15)) +
  annotate("text", x = min(mydata$X0) +2, y = min(mydata$X1) -1, label = "UMAP_1",
           color="black",size = 5 ) + 
  annotate("text", x = min(mydata$X0) -1, y = min(mydata$X1) + 2, label = "UMAP_2",
           color="black",size = 5 ,angle=90)


## 只含有免疫细胞的年轻人
young_data<-subset(mydata, level %in% 'adult')
young_data <- subset(young_data, cell_type2 %in% c("immune"))
nrow(young_data)
p3 <- ggplot2::ggplot(young_data,aes(x = X0, y = X1, color = cell_type))+
  geom_point(size = 0.1, alpha = 1) +
  scale_color_manual(values = c("#91D1C2FF","#B09C85FF","#8491B4FF","#DC0000FF","#00A087FF","#4DBBD5FF","#F39B7FFF","#3C5488FF"),
                     labels=c("TAM 2" = "microglia","TAM 1" = "MDM","prol. TAM" = "TAM 3"))+
  #scale_color_manual(values = c("brown3","chartreuse3","khaki1","plum1","lightskyblue","darkorange","hotpink","seashell3","#E87796","#E64B35FF"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+ #设置legend中 点的大小 
  ggtitle("Immune Cells in Adult") +
  cluster_theme +
  #绘制左下角箭头
  geom_segment(aes(x = min(mydata$X0) , y = min(mydata$X1) ,
                   xend = min(mydata$X0) +6, yend = min(mydata$X1) ),
               colour = "black", size=0.3,arrow = arrow(length = unit(0.3,"cm"),type="closed",angle=15))+ 
  geom_segment(aes(x = min(mydata$X0)  , y = min(mydata$X1)  ,
                   xend = min(mydata$X0) , yend = min(mydata$X1) + 6),
               colour = "black", size=0.3,arrow = arrow(length = unit(0.3,"cm"),type="closed",angle=15)) +
  annotate("text", x = min(mydata$X0) +2, y = min(mydata$X1) -1, label = "UMAP_1",
           color="black",size = 5 ) + 
  annotate("text", x = min(mydata$X0) -1, y = min(mydata$X1) + 2, label = "UMAP_2",
           color="black",size = 5 ,angle=90)

## 只含有免疫细胞的老年人
older_data<-subset(mydata, level %in% 'aged')
older_data <- subset(older_data, cell_type2 %in% c("immune"))
nrow(older_data)
p4 <- ggplot2::ggplot(older_data,aes(x = X0, y = X1, color = cell_type))+
  geom_point(size = 0.1, alpha = 1) +
  scale_color_manual(values = c("#91D1C2FF","#B09C85FF","#8491B4FF","#DC0000FF","#00A087FF","#4DBBD5FF","#F39B7FFF","#3C5488FF"),
                     labels=c("TAM 2" = "microglia","TAM 1" = "MDM","prol. TAM" = "TAM 3"))+
  #scale_color_manual(values = c("brown3","chartreuse3","khaki1","plum1","lightskyblue","darkorange","hotpink","seashell3","#E87796","#E64B35FF"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+ #设置legend中 点的大小 
  ggtitle("Immue Cells in Aged") +
  cluster_theme +
  #绘制左下角箭头
  geom_segment(aes(x = min(mydata$X0) , y = min(mydata$X1) ,
                   xend = min(mydata$X0) +6, yend = min(mydata$X1) ),
               colour = "black", size=0.3,arrow = arrow(length = unit(0.3,"cm"),type="closed",angle=15))+ 
  geom_segment(aes(x = min(mydata$X0)  , y = min(mydata$X1)  ,
                   xend = min(mydata$X0) , yend = min(mydata$X1) + 6),
               colour = "black", size=0.3,arrow = arrow(length = unit(0.3,"cm"),type="closed",angle=15)) +
  annotate("text", x = min(mydata$X0) +2, y = min(mydata$X1) -1, label = "UMAP_1",
           color="black",size = 5 ) + 
  annotate("text", x = min(mydata$X0) -1, y = min(mydata$X1) + 2, label = "UMAP_2",
           color="black",size = 5 ,angle=90)

setwd('/home/wsa/result/')
pdf('Cells in Adult.pdf',width = 8,height = 6)
p1
dev.off()
pdf('Cells in Aged.pdf',width = 8,height = 6)
p2
dev.off()
pdf('Immune Cells in Adult.pdf',width = 8,height = 6)
p3
dev.off()
pdf('Immune Cells in Aged.pdf',width = 8,height = 6)
p4
dev.off()

png("Cells in Adult.png",units="in", width=8, height=6,res=300)
p1
dev.off()
png("Cells in Aged.png",units="in", width=8, height=6,res=300)
p2
dev.off()
png("Immune Cells in Adult.png",units="in", width=8, height=6,res=300)
p3
dev.off()
png("Immune Cells in Aged.png",units="in", width=8, height=6,res=300)
p4
dev.off()

################################################################################
## 亚聚类
################################################################################
## 所有年龄层
table(subumap_data$level)
ggplot2::ggplot(subumap_data,aes(x = X0, y = X1, color = factor(sublabel)))+
  geom_point(size = 0.8, alpha = 1) +
  scale_color_manual(values = c("brown3","chartreuse3","khaki1","plum1","lightskyblue","darkorange","hotpink","seashell3","#E87796","#FE8980"))+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"),
        legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=20), #设置legend标签的大小
        legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5))) +  #设置legend中 点的大小 
  ggtitle("all") +
  theme(plot.title = element_text(hjust = 0.5)) #设置标题居中 #设置标题居中

## 未成年
paediatric_data<-subset(subumap_data, level %in% 'paediatric')
nrow(paediatric_data)
ggplot2::ggplot(paediatric_data,aes(x = X0, y = X1, color = factor(sublabel)))+
  geom_point(size = 0.8, alpha = 1) +
  scale_color_manual(values = c("brown3","chartreuse3","khaki1","plum1","lightskyblue","darkorange","hotpink","seashell3","#E87796","#FE8980"))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))+
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=20), #设置legend标签的大小
    legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5))) +  #设置legend中 点的大小 
  ggtitle("paediatric") +
  theme(plot.title = element_text(hjust = 0.5)) #设置标题居中

## 年轻人
young_data<-subset(subumap_data, level %in% 'young')
nrow(young_data)
ggplot2::ggplot(young_data,aes(x = X0, y = X1, color = factor(sublabel)))+
  geom_point(size = 0.8, alpha = 1) +
  scale_color_manual(values = c("brown3","chartreuse3","khaki1","plum1","lightskyblue","darkorange","hotpink","seashell3","#E87796","#FE8980"))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))+
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=20), #设置legend标签的大小
    legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5))) +  #设置legend中 点的大小 
  ggtitle("young") +
  theme(plot.title = element_text(hjust = 0.5)) #设置标题居中


## 老年人
older_data<-subset(subumap_data, level %in% 'older')
nrow(older_data)
ggplot2::ggplot(older_data,aes(x = X0, y = X1, color = factor(sublabel)))+
  geom_point(size = 0.8, alpha = 1) +
  scale_color_manual(values = c("brown3","chartreuse3","khaki1","plum1","lightskyblue","darkorange","hotpink","seashell3","#E87796","#FE8980"))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))+
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=20), #设置legend标签的大小
    legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5))) +  #设置legend中 点的大小 
  ggtitle("older") +
  theme(plot.title = element_text(hjust = 0.5)) #设置标题居中





################################################################################
## 细胞比例图
################################################################################
## 计算不同年龄每种免疫细胞比例,免疫细胞分组条形图
table(mydata$level)#查看各组细胞数
prop.table(table(mydata$cell_type))
table(mydata$cell_type, mydata$level)#各组不同细胞群细胞数
number <- table(mydata$cell_type, mydata$level)[c(1,2,3,4,5,6,8,9),c(1,2)]
number
Cellratio <- prop.table(number, margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)
Cellratio
colourCount = length(unique(Cellratio$Var1))

display.brewer.all()
display.brewer.pal(11,"Spectral")
color1 <- brewer.pal(11,"Spectral")
mycolor2 <- c(color1[9],color1[4],color1[3])

mycolor2 <- c("#606c38","#dda15e","#283618")
mycolor2 <- c("#264653","#2a9d8f","#f4a261")
mycolor2 <- c("#293241","#ee6c4d")
p1 <- ggplot(data = Cellratio,aes(x = Var1,y = Freq,fill = Var2))+
  geom_col(position = 'dodge',
           width = 0.5)+
  theme(panel.background = element_blank(),  #去除灰色背景
        panel.grid = element_blank(), #去除网格线
        axis.line = element_line(colour = "#000000",size = 0.3),
        axis.text = element_text(colour = "#000000" ,size = 18),
        axis.text.x = element_text(angle = 60,vjust = 1,hjust = 1), 
        axis.ticks = element_line(colour = "#000000" ,size = 0.3),
        axis.ticks.length = unit(1,'mm'),
        plot.margin = unit(c(0.5,0.4,0.4,0.3),"cm"),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 24,hjust = 0.5),
        legend.text = element_text(size=16), #设置legend标签的大小
        legend.key.size=unit(0.8,'cm')  ,  # 设置legend标签之间的大小
  )+
  scale_fill_manual(values=rev(mycolor2),
                    labels=c("Adult","Aged")) +
  scale_x_discrete(labels = c("TAM 1" = "MDM","TAM 2" = "microglia", "prol. TAM" = "TAM 3")) +
  #scale_fill_manual(values=c("#FBD84A","#666666")) +     #设置填充颜色
  labs(x=' ',y = 'Ratio')+  #设置横纵坐标名字
  guides(fill = guide_legend(title = NULL))+  # 删掉图例名称
  scale_y_continuous(labels=percent) +   # 使纵坐标呈现百分比
  ggtitle("Ratio of Cells")
setwd('/home/wsa/result/')
pdf('Ratio of Immune Cells.pdf',width = 8,height = 6)
p1
dev.off()
png("Ratio of Immune Cells.png",units="in", width=8, height=6,res=300)
p1
dev.off()


## 肿瘤占总细胞的条形图
subdata <- subset(mydata,!(group %in% '4'))
subdata <- subset(subdata,!(level %in% c('NA')))
table(subdata$level)#查看各组细胞数
prop.table(table(subdata$cell_type2))
table(subdata$cell_type2, subdata$level)#各组不同细胞群细胞数
number <- table(subdata$cell_type2, subdata$level)[,c(1,2)]
number
Cellratio <- prop.table(number, margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)
Cellratio

mycolor2 <- c("#264653","#f4a261")
mycolor2 <- c("#3C5488FF","#E64B35FF")
p1 <- ggplot(data = Cellratio,aes(x = Var1,y = Freq,fill = Var2))+
  geom_col(position = 'dodge',
           width = 0.5)+
  theme(panel.background = element_blank(),  #去除灰色背景
        panel.grid = element_blank(), #去除网格线
        axis.line = element_line(colour = "#000000",size = 0.3),
        axis.text = element_text(colour = "#000000" ,size = 18),
        axis.text.x = element_text(angle = 60,vjust = 1,hjust = 1), 
        axis.ticks = element_line(colour = "#000000" ,size = 0.3),
        axis.ticks.length = unit(1,'mm'),
        plot.margin = unit(c(0.5,0.4,0.4,0.3),"cm"),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 24,hjust = 0.5),
        legend.text = element_text(size=16), #设置legend标签的大小
        legend.key.size=unit(0.8,'cm')  ,  # 设置legend标签之间的大小
  )+
  scale_fill_manual(values=rev(mycolor2),
                    labels=c("Adult","Aged")) +
  scale_x_discrete(labels = c("TAM 1" = "MDM","TAM 2" = "microglia", "prol. TAM" = "TAM 3")) +
  #scale_fill_manual(values=c("#FBD84A","#666666")) +     #设置填充颜色
  labs(x=' ',y = 'Ratio')+  #设置横纵坐标名字
  guides(fill = guide_legend(title = NULL))+  # 删掉图例名称
  scale_y_continuous(labels=percent) +   # 使纵坐标呈现百分比
  ggtitle("Ratio of Cells")
setwd('/home/wsa/result/')
pdf('Ratio of Cells.pdf',width = 8,height = 6)
p1
dev.off()
png("Ratio of Cells.png",units="in", width=8, height=6,res=300)
p1
dev.off()



## microglia占总细胞的条形图
group <- seurat@meta.data$group
mydata$group <- group
subdata <- subset(mydata,!(group %in% '4'))
table(subdata$level)#查看各组细胞数
prop.table(table(subdata$cell_type))
table(subdata$cell_type, subdata$level)#各组不同细胞群细胞数
number <- table(subdata$cell_type, subdata$level)[,c(1,2)]
number
Cellratio <- prop.table(number, margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)
Cellratio

mycolor2 <- c("#293241","#ee6c4d","#98c1d9")
ggplot(data = Cellratio,aes(x = Var1,y = Freq,fill = Var2))+
  geom_col(position = 'dodge',
           width = 0.5)+
  theme(panel.background = element_blank(),  #去除灰色背景
        panel.grid = element_blank(), #去除网格线
        axis.line = element_line(colour = "#000000",size = 0.3),
        axis.text = element_text(colour = "#000000" ,size = 13),
        axis.text.x = element_text(angle = 60,vjust = 1,hjust = 1), 
        axis.ticks = element_line(colour = "#000000" ,size = 0.3),
        axis.ticks.length = unit(1,'mm'),
        plot.margin = unit(c(0.5,0.4,0.4,0.3),"cm"),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 20,hjust = 0.5),
  )+
  scale_fill_manual(values=rev(mycolor2),
                    labels=c("Adult","Aged")) +
  scale_x_discrete(labels = c("TAM 1" = "MDM","TAM 2" = "microglia", "prol. TAM" = "TAM 3")) +
  #scale_fill_manual(values=c("#FBD84A","#666666")) +     #设置填充颜色
  labs(x=' ',y = 'Ratio')+  #设置横纵坐标名字
  guides(fill = guide_legend(title = NULL))+  # 删掉图例名称
  scale_y_continuous(labels=percent) +   # 使纵坐标呈现百分比
  ggtitle("ratio of all cells")



## TAM1等亚群分组条形图
table(subumap_data$level)#查看各组细胞数
prop.table(table(subumap_data$sublabel))
table(subumap_data$sublabel, subumap_data$level)#各组不同细胞群细胞数
number <- table(subumap_data$sublabel, subumap_data$level)[,c(2,3,1)]
number
Cellratio <- prop.table(number, margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)
Cellratio
colourCount = length(unique(Cellratio$Var1))

mycolor2 <- c("#293241","#ee6c4d","#98c1d9")
ggplot(data = Cellratio,aes(x = Var1,y = Freq,fill = Var2))+
  geom_col(position = 'dodge',
           width = 0.7)+
  theme(panel.background = element_blank(),  #去除灰色背景
        panel.grid = element_blank(), #去除网格线
        axis.line = element_line(colour = "#000000",size = 0.3),
        axis.text = element_text(colour = "#000000" ,size = 13),
        axis.text.x = element_text(angle = 60,vjust = 1,hjust = 1), 
        axis.ticks = element_line(colour = "#000000" ,size = 0.3),
        axis.ticks.length = unit(1,'mm'),
        plot.margin = unit(c(0.5,0.4,0.4,0.3),"cm"),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 20,hjust = 0.5),
  )+
  scale_fill_manual(values=rev(mycolor2),
                    name="Experimental",
                    breaks=c("paediatric","young","older"),
                    labels=c("child","Adult","Aged")) +
  #scale_fill_manual(values=c("#FBD84A","#666666")) +     #设置填充颜色
  labs(x=' ',y = 'Ratio')+  #设置横纵坐标名字
  guides(fill = guide_legend(title = NULL))+  # 删掉图例名称
  scale_y_continuous(labels=percent) +   # 使纵坐标呈现百分比
  ggtitle("ratio of  subcluster in TAM3")



##绘制不同年龄堆叠柱状图
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))


################################################################################
## 散点图
################################################################################
### 免疫细胞占总细胞
subdata <- subset(mydata,!(group %in% '4'))
table(subdata$sample)#查看各组细胞数
prop.table(table(subdata$cell_type2))
table(subdata$cell_type2, subdata$sample)#各组不同细胞群细胞数
number <- table(subdata$cell_type2, subdata$sample)
number
Cellratio <- prop.table(number, margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- Cellratio[1,]
Cellratio
Cellratio <- as.data.frame(Cellratio)
Cellratio
txt <- read.table('/home/wsa/data/plot/sample_age_1850_6590.txt',sep='\t',header = TRUE)
sample = txt$sample
for (j in sample){
  #anno_rownumber <- which(rownames(Cellratio) %in% j, arr.ind = TRUE)
  txt_rownumber <- which(txt$sample %in% j, arr.ind = TRUE)
  Cellratio[j,'age'] <- txt[txt_rownumber,"age"]
  Cellratio[j,'level'] <- txt[txt_rownumber,"level"]
}
data = na.omit(Cellratio)
ggplot(data=data, # 定义源数据
       mapping=aes(
                x=age,
                y=Cellratio))+
  geom_point()+  # 以散点图形式呈现
  geom_smooth()   #拟合曲线


### 小胶质细胞占免疫细胞
txt <- read.table('/home/wsa/data/plot/sample_age_1850_6590.txt',sep='\t',header = TRUE)
table(mydata$sample)#查看各组细胞数
prop.table(table(mydata$cell_type))
table(mydata$cell_type, mydata$sample)#各组不同细胞群细胞数
number <- table(mydata$cell_type, mydata$sample)[c(1,2,3,4,6,7,8,9),]
number
Cellratio <- prop.table(number, margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- Cellratio[8,]
Cellratio
Cellratio <- as.data.frame(Cellratio)
Cellratio
sample = txt$sample
for (j in sample){
  txt_rownumber <- which(txt$sample %in% j, arr.ind = TRUE)
  Cellratio[j,'age'] <- txt[txt_rownumber,"age"]
  Cellratio[j,'level'] <- txt[txt_rownumber,"level"]
}
data = na.omit(Cellratio)
ggplot(data=data, # 定义源数据
       mapping=aes(
         x=age,
         y=Cellratio))+
  geom_point()+  # 以散点图形式呈现
  geom_smooth()   #拟合曲线



### MDM占免疫细胞
txt <- read.table('/home/wsa/data/plot/sample_age_1850_6590.txt',sep='\t',header = TRUE)
table(mydata$sample)#查看各组细胞数
prop.table(table(mydata$cell_type))
table(mydata$cell_type, mydata$sample)#各组不同细胞群细胞数
number <- table(mydata$cell_type, mydata$sample)[c(1,2,3,4,6,7,8,9),]
number
Cellratio <- prop.table(number, margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- Cellratio[7,]
Cellratio
Cellratio <- as.data.frame(Cellratio)
Cellratio
sample = txt$sample
for (j in sample){
  txt_rownumber <- which(txt$sample %in% j, arr.ind = TRUE)
  Cellratio[j,'age'] <- txt[txt_rownumber,"age"]
  Cellratio[j,'level'] <- txt[txt_rownumber,"level"]
}
data = na.omit(Cellratio)
ggplot(data=data, # 定义源数据
       mapping=aes(
         x=age,
         y=Cellratio))+
  geom_point()+  # 以散点图形式呈现
  geom_smooth()   #拟合曲线



################################################################################
## microglia/immune
################################################################################
txt <- read.table('/home/wsa/data/plot/sample_age_1850_6590.txt',sep='\t',header = TRUE)
table(mydata$sample)#查看各组细胞数
prop.table(table(mydata$cell_type))
table(mydata$cell_type, mydata$sample)#各组不同细胞群细胞数
number <- table(mydata$cell_type, mydata$sample)[c(1,2,3,4,5,6,8,9),]
number
Cellratio <- prop.table(number, margin = 2)#计算各组样本不同细胞群比例
Cellratio
setwd('/home/wsa/data/cell_ratio')
write.table(Cellratio, file = 'only_immune.txt',sep='\t')
Cellratio <- as.data.frame(Cellratio)
Cellratio
sample = txt$sample
for (j in sample){
  anno_rownumber <- which(Cellratio$Var2 %in% j, arr.ind = TRUE)
  txt_rownumber <- which(txt$sample %in% j, arr.ind = TRUE)
  Cellratio[anno_rownumber,'age'] <- txt[txt_rownumber,"age"]
  Cellratio[anno_rownumber,'level'] <- txt[txt_rownumber,"level"]
}
order(Cellratio$age)
Cellratio <- Cellratio[order(Cellratio$age),]
View(Cellratio)
data = na.omit(Cellratio)

faceted_label <- c(
  `18-50` = 'Adult',
  `65-100` = 'Aged')
display.brewer.all()
color1 <- brewer.pal(11,"Spectral")
color2 <- brewer.pal(11,"RdGy")
color3 <- brewer.pal(11,"BrBG")
mycolor2 <- c(color3[3],color1[11],color1[3],color1[10],color2[11],color2[10],color2[9],color2[8])
mycolor2 <- c("#264653","#2a9d8f","#e9c46a","#f4a261","#e76f51",color2[10],color2[9],color2[8])
mycolor2 <- c("#606c38","#283618","#fefae0","#dda15e","#bc6c25",color2[10],color2[9],color2[8])
mycolor2 <- c("#293241","#ee6c4d","#e0fbfc","#98c1d9","#3d5a80",color2[10],color2[9],color2[8])
p <- ggplot(data) + 
  geom_bar(aes(x =reorder(Var2,age), y= Freq, fill = Var1),stat = "identity",width = 0.8,size = 0.5,position = 'fill')+  #reorder函数，让x轴按照年龄顺序排布
  theme(panel.background = element_blank(),  #去除灰色背景
        panel.grid = element_blank(), #去除网格线
        axis.line = element_line(colour = "#000000",size = 0.3),
        axis.text = element_text(colour = "#000000" ),
        axis.text.x = element_text(size = 18, angle = 60,vjust = 1,hjust = 1), 
        axis.text.y = element_text(size = 18),
        axis.ticks = element_line(colour = "#000000" ,size = 0.3),
        axis.ticks.length = unit(1,'mm'),
        plot.margin = unit(c(0.5,0.4,0.4,0.3),"cm"),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 24,hjust = 0.5),
        legend.text = element_text(size=16), #设置legend标签的大小
        legend.key.size=unit(0.8,'cm')  ,  # 设置legend标签之间的大小
        strip.background.x = element_rect(fill = "#303030", colour = "#303030"),
        strip.text.x = element_text(size = 20, face="bold",colour = "#FFFFFF")
  )+ 
  scale_fill_manual(values=rev(mycolor2),
                    labels=c("TAM 2" = "microglia","TAM 1" = "MDM","prol. TAM" = "TAM 3"))+
  #scale_fill_manual(values=c("#FBD84A","#666666")) +     #设置填充颜色
  labs(x=' ',y = 'Ratio')+  #设置横纵坐标名字
  guides(fill = guide_legend(title = NULL))+  # 删掉图例名称
  scale_y_continuous(labels=percent) +   # 使纵坐标呈现百分比
  ggtitle("Ratio of Immune Cells in every samples") +
  facet_grid(.~level,scales='free_x',space='free_x',labeller = as_labeller(faceted_label))
setwd('/home/wsa/result/')
pdf('every samples.pdf',width = 8,height = 6)
p
dev.off()
png("every samples.png",units="in", width=8, height=6,res=300)
p
dev.off()


################################################################################
## tumor/total 每个样本肿瘤细胞比例,分组柱状图
################################################################################
subdata <- subset(mydata,!(group %in% c('4')))
subdata <- subset(subdata,!(level %in% c('NA')))
table(subdata$sample)#查看各组细胞数
table(subdata$level)
prop.table(table(subdata$cell_type2))
table(subdata$cell_type2, subdata$sample)#各组不同细胞群细胞数
number <- table(subdata$cell_type2, subdata$sample)
number
Cellratio <- prop.table(number, margin = 2)#计算各组样本不同细胞群比例
Cellratio
setwd('/home/wsa/data/cell_ratio')
write.table(Cellratio, file = 'patient.txt',sep='\t')





################################################################################
## microglia/total 每个样本肿瘤细胞比例
################################################################################
subdata <- subset(mydata,!(group %in% c('4')))
subdata <- subset(subdata,!(level %in% c('NA')))
table(subdata$level)
table(subdata$sample)#查看各组细胞数
prop.table(table(subdata$cell_type))
table(subdata$cell_type, subdata$sample)#各组不同细胞群细胞数
number <- table(subdata$cell_type, subdata$sample)
number
Cellratio <- prop.table(number, margin = 2)#计算各组样本不同细胞群比例
Cellratio
setwd('/home/wsa/data/cell_ratio')
write.table(Cellratio, file = 'microglia_total.txt',sep='\t')



Cellratio <- as.data.frame(Cellratio)
Cellratio
txt <- read.table('/home/wsa/data/plot/sample_age_1850_6590.txt',sep='\t',header = TRUE)
sample = txt$sample
for (j in sample){
  anno_rownumber <- which(Cellratio$Var2 %in% j, arr.ind = TRUE)
  txt_rownumber <- which(txt$sample %in% j, arr.ind = TRUE)
  Cellratio[anno_rownumber,'age'] <- txt[txt_rownumber,"age"]
  Cellratio[anno_rownumber,'level'] <- txt[txt_rownumber,"level"]
}
order(Cellratio$age)
Cellratio <- Cellratio[order(Cellratio$age),]
View(Cellratio)
data = na.omit(Cellratio)

faceted_label <- c(
  `18-50` = 'Adult',
  `65-100` = 'Aged')
mycolor2 <- c("#606c38","#283618","#fefae0","#dda15e","#bc6c25",color2[10],color2[9],color2[8])
mycolor2 <- c("#293241","#ee6c4d","#e0fbfc","#98c1d9","#3d5a80","#dda15e","#bc6c25",color2[10],color2[9],color2[8])
p <- ggplot(data) + 
  geom_bar(aes(x =reorder(Var2,age), y= Freq, fill = Var1),stat = "identity",width = 0.8,size = 0.5,position = 'fill')+  #reorder函数，让x轴按照年龄顺序排布
  theme(panel.background = element_blank(),  #去除灰色背景
        panel.grid = element_blank(), #去除网格线
        axis.line = element_line(colour = "#000000",size = 0.3),
        axis.text = element_text(colour = "#000000" ),
        axis.text.x = element_text(size = 8, angle = 60,vjust = 1,hjust = 1), 
        axis.text.y = element_text(size = 12),
        axis.ticks = element_line(colour = "#000000" ,size = 0.3),
        axis.ticks.length = unit(1,'mm'),
        plot.margin = unit(c(0.5,0.4,0.4,0.3),"cm"),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 20,hjust = 0.5),
        strip.background.x = element_rect(fill = "#303030", colour = "#303030"),
        strip.text.x = element_text(size = 10, face="bold",colour = "#FFFFFF")
  )+ 
  scale_fill_manual(values=rev(mycolor2),
                    labels=c("TAM 2" = "microglia","TAM 1" = "MDM","prol. TAM" = "TAM 3"))+
  #scale_fill_manual(values=c("#FBD84A","#666666")) +     #设置填充颜色
  labs(x=' ',y = 'Ratio')+  #设置横纵坐标名字
  guides(fill = guide_legend(title = NULL))+  # 删掉图例名称
  scale_y_continuous(labels=percent) +   # 使纵坐标呈现百分比
  ggtitle("ratio of all cells in every sample") +
  facet_grid(.~level,scales='free_x',space='free_x',labeller = as_labeller(faceted_label))
p

if (FALSE) {
##每个样本肿瘤细胞比例,按年龄大小排好顺序，去除group4的样本（没有肿瘤）
mydata$cell_type2 <- mydata$cell_type
mydata[which(!(mydata$cell_type2 %in% c("NormalBrain","Tumor"))),'cell_type2'] <- 'immune'
subdata <- subset(mydata,!(sample %in% c('ND1','ND2','ND4','ND5','ND7')))
table(subdata$sample)#查看各组细胞数
prop.table(table(subdata$cell_type2))
table(subdata$cell_type2, subdata$sample)#各组不同细胞群细胞数
number <- table(subdata$cell_type2,subdata$sample)
number
Cellratio <- prop.table(number, margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)
Cellratio
child <- c('BT749',
       'BT771',
       'BT830',
       'BT920',
       'BT1160',
       'MUV1',
       'MUV5',
       'MUV10',
       'BCH836',
       'BCH869',
       'BCH1126')  
anno_rownumber <- which(Cellratio$Var2 %in% child, arr.ind = TRUE)
Cellratio[anno_rownumber,'level'] <- '0-18'
adult <- c('S1',
           'S2',
           'S11',
           'ND7',
           '8036')
anno_rownumber <- which(Cellratio$Var2 %in% adult, arr.ind = TRUE)
Cellratio[anno_rownumber,'level'] <- '30-50'
aged <- c('S3',
          'S5',
          'S12',
          'MGH102',
          'MGH105',
          'MGH115',
          'MGH124',
          'MGH126',
          'MGH143',
          'MGH102',
          'MGH104',
          'MGH105',
          'MGH106',
          'MGH110',
          'MGH113',
          'MGH115',
          'MGH121',
          'MGH122',
          'MGH124',
          'MGH128',
          'MGH143',
          'MGH151',
          'ND1',
          'ND2',
          'ND4',
          'ND5',
          '8165C',
          '8165PV')
anno_rownumber <- which(Cellratio$Var2 %in% aged, arr.ind = TRUE)
Cellratio[anno_rownumber,'level'] <- '60-70'
data = na.omit(Cellratio)

faceted_label <- c(
  `0-18` = 'Child',
  `30-50` = 'Adult',
  `60-70` = 'Aged')
mycolor2 <- c("#293241","#ee6c4d","#98c1d9")
data$Var1 = factor(data$Var1, levels = c("NormalBrain","immune","Tumor")) #调整填充顺序
ggplot(data) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.8,size = 0.5,position = 'fill')+  #reorder函数，让x轴按照年龄顺序排布
  theme(panel.background = element_blank(),  #去除灰色背景
        panel.grid = element_blank(), #去除网格线
        axis.line = element_line(colour = "#000000",size = 0.3),
        axis.text = element_text(colour = "#000000" ),
        axis.text.x = element_text(size = 8, angle = 60,vjust = 1,hjust = 1), 
        axis.text.y = element_text(size = 12),
        axis.ticks = element_line(colour = "#000000" ,size = 0.3),
        axis.ticks.length = unit(1,'mm'),
        plot.margin = unit(c(0.5,0,0.4,0.4),"cm"),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 20,hjust = 0.5),
        strip.background.x = element_rect(fill = "#303030", colour = "#303030"),
        strip.text.x = element_text(size = 10, face="bold",colour = "#FFFFFF")
  )+ 
  scale_fill_manual(values=rev(mycolor2),
                    name="",
                    limits = c("NormalBrain","immune","Tumor"))+
  #scale_fill_manual(values=c("#FBD84A","#666666")) +     #设置填充颜色
  labs(x=' ',y = 'Ratio')+  #设置横纵坐标名字
  guides(fill = guide_legend(title = NULL))+  # 删掉图例名称
  scale_y_continuous(labels=percent) +   # 使纵坐标呈现百分比
  ggtitle("ratio of cells in every sample") +
  facet_grid(.~level,scales='free_x',space='free_x',labeller = as_labeller(faceted_label))
}




## 计算未成年人每个样本免疫细胞比例
table(paediatric_data$sample)#查看各组细胞数
prop.table(table(paediatric_data$cell_type))
table(paediatric_data$cell_type, paediatric_data$sample)#各组不同细胞群细胞数
number <- table(paediatric_data$cell_type, paediatric_data$sample)[c(1,2,3,4,6,7,8,9),]
number
Cellratio <- prop.table(number, margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)
Cellratio
colourCount = length(unique(Cellratio$Var1))
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.8,size = 0.5,position = 'fill')+ 
  theme(panel.background = element_blank(),  #去除灰色背景
      panel.grid = element_blank(), #去除网格线
      axis.line = element_line(colour = "#000000",size = 0.3),
      axis.text = element_text(colour = "#000000" ,size = 13),
      axis.text.x = element_text(angle = 40,vjust = 0.5,hjust = 0.5), 
      axis.ticks = element_line(colour = "#000000" ,size = 0.3),
      axis.ticks.length = unit(1,'mm'),
      plot.margin = unit(c(0.5,0,0,0),"cm"),
      axis.title.y = element_text(size = 15),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 20,hjust = 0.5),
)+
  #scale_fill_manual(values=c("#FBD84A","#666666")) +     #设置填充颜色
  labs(x=' ',y = 'Ratio')+  #设置横纵坐标名字
  guides(fill = guide_legend(title = NULL))+  # 删掉图例名称
  scale_y_continuous(labels=percent) +   # 使纵坐标呈现百分比
  ggtitle("ratio of immune cells in every sample") +
  theme(plot.title = element_text(hjust = 0.5))   #设置标题居中





## 计算年轻人每个样本免疫细胞比例
table(young_data$sample)#查看各组细胞数
prop.table(table(young_data$cell_type))
table(young_data$cell_type, young_data$sample)#各组不同细胞群细胞数
number <- table(young_data$cell_type, young_data$sample)[c(1,2,3,4,6,7,8,9),]
number
Cellratio <- prop.table(number, margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)
Cellratio
colourCount = length(unique(Cellratio$Var1))
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  ggtitle("young") +
  theme(plot.title = element_text(hjust = 0.5)) #设置标题居中

## 计算老年人每个样本免疫细胞比例
table(older_data$sample)#查看各组细胞数
prop.table(table(older_data$cell_type))
table(older_data$cell_type, older_data$sample)#各组不同细胞群细胞数
number <- table(older_data$cell_type, older_data$sample)[c(1,2,3,4,6,7,8,9),]
number
Cellratio <- prop.table(number, margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)
Cellratio
colourCount = length(unique(Cellratio$Var1))
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  ggtitle("older") +
  theme(plot.title = element_text(hjust = 0.5)) #设置标题居中




##绘制饼图
##http://www.sthda.com/english/wiki/ggplot2-pie-chart-quick-start-guide-r-software-and-data-visualization
young_cellratio <- subset(Cellratio,subset = (Cellratio$Var2=='young'))
older_cellratio <- subset(Cellratio,subset = (Cellratio$Var2=='older'))
#年轻人
# Barplot
bp<- ggplot(young_cellratio, aes(x="", y=Freq, fill=Var1))+
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0)
#pie + scale_fill_manual(values=c("#2F4858","#999999","#106397","#5C69A9", "#946CAE","#C46FA7","#E87796","#FE8980"))
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
# Apply blank theme
library(scales)
pie +  blank_theme +
  theme(axis.text.x=element_blank())+
  ggtitle("young") +
  theme(plot.title = element_text(hjust = 0.5)) #设置标题居中

#老年人
# Barplot
bp<- ggplot(older_cellratio, aes(x="", y=Freq, fill=Var1))+
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0)
#pie + scale_fill_manual(values=c("#2F4858","#999999","#106397","#5C69A9", "#946CAE","#C46FA7","#E87796","#FE8980"))
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
# Apply blank theme
library(scales)
pie +  blank_theme +
  theme(axis.text.x=element_blank())+
  ggtitle("older") +
  theme(plot.title = element_text(hjust = 0.5)) #设置标题居中