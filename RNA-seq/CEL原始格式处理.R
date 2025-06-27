# 本视频以GEO数据库GSE107731数据集为例
library(affyPLM)
library(limma)
##### 设置CEL数据路径 ####
# 根据个人数据设置工作目录
setwd("/data/gpfs02/wmo/work/HF/RNA_seq/LW/rawdata")

##### 读取CEL数据 ####
Data <- ReadAffy()

##### 预处理的一体化算法 ####
# MAS5和RMA算法使用较多；
# RMA处理后的数据是经过2为底的对数转换的，而MAS5不是；
# 大多数芯片分析软件和函数的输入数据必须经过对数变换。
# eset.rma <- rma(Data)
eset.mas5 <- mas5(Data)
eset.mas5 <- log2(eset.mas5 + 1)

##### 提取一体化算法计算表达量 ####
exp.rma <- exprs(eset.rma)

##### 预处理表达矩阵 ####
exp.rma <- cbind(rownames(exp.rma), exp.rma)
exp.rma <- as.data.frame(exp.rma)
colnames(exp.rma)[2:ncol(exp.rma)] <- substr(colnames(exp.rma)[2:ncol(exp.rma)], 1, 10)

##### probid转换 ####
setwd("/data/gpfs02/wmo/work/HF/RNA_seq/LW")
plt <- read.table("GPL1261-56135.txt", sep = "\t", header = T)
exp.rma <- exp.rma[exp.rma[, 1] %in% plt[, 1], ]
exp.rma[, 1] <- plt[which(exp.rma[, 1] %in% plt[, 1]), 2]

##### 相同基因表达取平均 ####
exp.rma <- as.matrix(exp.rma) # 转换为矩阵
rownames(exp.rma) <- exp.rma[, 1] # 将第一列设置为行名
exp.rma <- exp.rma[, 2:ncol(exp.rma)] # 剔除第一行行名
exp.rma <- avereps(exp.rma) # 表达矩阵内重复基因的值替换为其平均值
exp.rma <- as.data.frame(exp.rma)
for (i in 1:ncol(exp.rma)) {
  exp.rma[, i] <- as.numeric(exp.rma[, i])
}

##### 读取分组数据 ####
group <- read.table("Group.txt", sep = "\t", header = T)

##### 分组矩阵 ####
class <- as.character(group$Group)
design <- model.matrix(~ 0 + factor(class))
colnames(design) <- c("disease", "normal")

##### 差异比较矩阵 ####
cont.matrix <- makeContrasts(normal - disease, levels = design)

##### 差异分析第一步 ####
fit <- lmFit(exp.rma, design)
fit <- contrasts.fit(fit, cont.matrix)

##### 差异分析第二步 ####
fit <- eBayes(fit)

##### 差异分析第三步 ####
analysis_result <- topTable(fit, adjust = "BH", number = 50000)

##### 提取差异表达基因 #####
diff_result <- analysis_result[(abs(analysis_result$logFC) >= 1 & analysis_result$adj.P.Val < 0.05), ]

##### 保存数据 ####
# write.table(analysis_result, "limma_analysis_result.txt", quote = F, sep = "\t")
write.table(diff_result, "limma_analysis_diff_result.txt", quote = F, sep = "\t")
