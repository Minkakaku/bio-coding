library(edgeR)
exp <- read.delim("all_exp.txt", row.names = 1)
exp.rma <- as.matrix(exp.rma)
exp.rma <- avereps(exp.rma) # 表达矩阵内重复基因的值替换为其平均值(可选)
head(exp)

d0 <- DGEList(exp)
# 注意： calcNormFactors并不会标准化数据，只是计算标准化因子
d0 <- calcNormFactors(d0)
d0

# 过滤低表达
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left

# sample names
snames <- colnames(exp)
snames
# 此数据有两个因素：cultivar（C,I5/I8）和time(6,9)
cultivar <- substr(snames, 1, nchar(snames) - 2) 
time <- substr(snames, nchar(snames) - 1, nchar(snames) - 1)
cultivar
time

# Create a new variable “group” that combines cultivar and time
group <- interaction(cultivar, time)
group

# Multidimensional scaling (MDS) plot
plotMDS(d, col = as.numeric(group))

mm <- model.matrix(~0 + factor(group))
colnames(mm) <- c("disease", "normal")

# plot Voom
par(mfrow = c(1, 3))
vd <- voom(d, mm, plot = T)
#voom的曲线应该很光滑，比较一下过滤低表达gene之前的图形：
voom(d0, mm, plot = T)
# lmFit fits a linear model using weighted least squares for each gene:
fit <- lmFit(vd, mm)
head(coef(fit))



cont.matrix <- makeContrasts(normal - disease, levels = colnames(coef(fit)))

fit <- contrasts.fit(fit, cont.matrix)

fit <- eBayes(fit)
plotSA(fit, main="Final model: Mean-variance trend")

analysis_result <- topTable(fit, adjust = "BH", number = 50000)

#提取差异表达基因
diff_result <- analysis_result[(abs(analysis_result$logFC) >= 1 & analysis_result$adj.P.Val < 0.05), ]

write.table(analysis_result, "limma_analysis_result.txt", quote = F, sep = "\t")
write.table(diff_result, "limma_analysis_diff_result.txt", quote = F, sep = "\t")
