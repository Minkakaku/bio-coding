rm(list=ls())
# 显示当前工作目录
cat("当前工作目录: ", getwd(), "\n")

# 设置工作目录，可根据需要修改路径。"." 表示当前目录，Windows 下使用正斜杠 /
workingDir = "."
setwd(workingDir)

# 加载 WGCNA 包
library(WGCNA)

# 重要设置，请勿省略
options(stringsAsFactors = FALSE)

# 启用 WGCNA 多线程，加速部分计算
# 注意：若使用 RStudio 或其他第三方 R 环境，请跳过此命令
# 当前此命令为代码运行所必需，若出现错误可忽略，但建议更新 WGCNA 包
if (!Sys.getenv("RSTUDIO") == "1") {
  enableWGCNAThreads()
}

# 加载第一部分保存的数据
lnames <- load(file = "V:\\repo\\bio-coding\\WGCNA\\opt_inputdata\\opt_inputdata.RData")
# 输出加载的变量名
cat("已加载的变量名: ", paste(lnames, collapse = ", "), "\n")

# 检查 datExpr 是否加载成功
if (!"datExpr" %in% lnames) {
  stop("未成功加载 datExpr 数据，请检查 'opt_inputdata.RData' 文件")
}

# 选择一组软阈值幂次
powers <- c(1:20)
powers
# 调用网络拓扑分析函数
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# 打开绘图窗口并设置布局
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 <- 0.9

# 绘制无标度拓扑拟合指数随软阈值幂次的变化
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
# 绘制 R^2 阈值线
abline(h = 0.90, col = "red")

# 绘制平均连接度随软阈值幂次的变化
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")

# 根据结果自动选择合适的软阈值
optimal_power <- sft$fitIndices[which.max(-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]), 1]
cat("选择的最优软阈值幂次为: ", optimal_power, "\n")

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

net <- blockwiseModules(datExpr, power = optimal_power,
                        TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "femaleMouseTOM", 
                        verbose = 3)

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

# 打开绘图窗口
sizeGrWindow(12, 9)

# 将标签转换为颜色以便绘图
mergedColors <- labels2colors(net$colors)

# 绘制树状图并在下方显示模块颜色
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]

# 保存分析结果
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "FemaleLiver-02-networkConstruction-auto.RData")
