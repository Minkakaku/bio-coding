##### 相同基因表达取平均 ####
exp.rma <- as.matrix(exp.rma) # 转换为矩阵
rownames(exp.rma) <- exp.rma[, 1] # 将第一列设置为行名
exp.rma <- exp.rma[, 2:ncol(exp.rma)] # 剔除第一行行名
exp.rma <- avereps(exp.rma) # 表达矩阵内重复基因的值替换为其平均值
exp.rma <- as.data.frame(exp.rma)
for (i in seq_len(length(ncol(exp.rma)))) {
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
