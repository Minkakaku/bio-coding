library(WGCNA)
library(reshape2)
library(stringr)

directory = 'D:/Users/24432/Desktop/����/WGCNA/csx'
setwd(directory)
list.files()

exprMat <- "All expression BG"

options(stringsAsFactors = FALSE)
# �򿪶��߳�
enableWGCNAThreads()

# �ٷ��Ƽ� "signed" �� "signed hybrid"
type = "unsigned"
# corType: pearson or bicor
corType = "pearson"

corFnc = ifelse(corType=="pearson", cor, bicor)
# �Զ�Ԫ��������������״��Ϣ���������ʱ��
# �����������������ڼ���״̬ʱ���������������
maxPOutliers = ifelse(corType=="pearson",1,0.05)

# ������Ʒ��״�Ķ�Ԫ����ʱ������
#robustY = ifelse(corType=="pearson",T,F)

dataExpr <- read.csv("All expression BG.csv", header = T, stringsAsFactors = F)
dataExpr2 <- dataExpr[-which(dataExpr$Other.Gene.ID == ""| dataExpr$Other.Gene.ID == "1-Mar"| dataExpr$Other.Gene.ID == "2-Mar"),]
rownames(dataExpr2) <- dataExpr2$Other.Gene.ID

dataExpr3 <- as.data.frame(sapply(dataExpr2[,c(2:17)],as.numeric))
dataExpr3[is.na(dataExpr3)] = 0
#dataExpr3 = dataExpr2[complete.cases(dataExpr2),]
#dataExpr3 <- dataExpr3[-which(rowMeans(dataExpr3) == dataExpr3$P1BGwt1.FPKM|rowMeans(dataExpr3)<5),]

# ɸѡ��λ����ƫ��ǰ75%�Ļ�������MAD����0.01
m.mad <- apply(dataExpr3,1,mad)
dataExprVar <- dataExpr3[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]

dataExpr4 <- as.data.frame(t(dataExprVar))

## ���ȱʧֵ
gsg = goodSamplesGenes(dataExpr4, verbose = 3)

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr4)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr4)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr4 = dataExpr4[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr4)
nSamples = nrow(dataExpr4)

## �鿴�Ƿ�����Ⱥ��Ʒ
sampleTree = hclust(dist(dataExpr4), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

## ����ֵ��ɸѡԭ����ʹ����������������ޱ����������
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr4, powerVector=powers, 
                        networkType=type, verbose=5)

par(mfrow = c(1,2))
cex1 = 0.9
# ������Soft threshold (power)���������ޱ�������������������ֵԽ�ߣ�
# ����Խ�����ޱ������ (non-scale)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# ɸѡ��׼��R-square=0.85 �������
abline(h=0.85,col="red")

# Soft threshold��ƽ����ͨ��
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
power = sft$powerEstimate
power

# �����ȷʵ���������������仯����ģ�Ҳ����ʹ������ľ���powerֵ��
  if (is.na(power)){
    power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                   ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                          ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                                 ifelse(type == "unsigned", 6, 12))       
                   )
    )
   }
power
##һ�������繹����One-step network construction and module detection##

net = blockwiseModules(dataExpr4, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0(exprMat, ".tom"),
                       verbose = 3)

# ����ģ���л�����Ŀ�Ķ��٣��������У����α��Ϊ `1-���ģ����`��
# **0 (grey)**��ʾ**δ**�����κ�ģ��Ļ��� 
table(net$colors)

## ��ɫ��Ϊ**δ����**��ģ��Ļ���
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# Plot the dendrogram and the module colors underneath
# ����Խ�������⣬������recutBlockwiseTrees����ʡ����ʱ��
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# module eigengene, ���Ի�����ͼ����Ϊÿ��ģ��Ļ���������Ƶ�չʾ
MEs = net$MEs

### ����Ҫ���¼��㣬���������־ͺ�
### �ٷ��̳������¼���ģ���ʼ���Բ�����ô�鷳
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

# ���ݻ������������о������õ��ĸ�ģ���������ͼ
# marDendro/marHeatmap �����¡����ϡ��ҵı߾�
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)

# ������÷ֲ����㣬�����õ�blocksize>=�ܻ�������ֱ��load����õ�TOM���
# ������Ҫ�ټ���һ�飬�ȽϺķ�ʱ��
# TOM = TOMsimilarityFromExpr(dataExpr, power=power, corType=corType, networkType=type)
load(net$TOMFiles[1], verbose=T)

## Loading objects:
##   TOM
#load("All expression BG.tom-block.1.RData")
TOM <- as.matrix(TOM)

dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong 
# connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function

# ��һ�����ر��ʱ������ͬʱ���㼶����
TOMplot(plotTOM, net$dendrograms, moduleColors, 
        main = "Network heatmap plot, all genes")
















