if (TRUE) {
    library(edgeR)
library(dplyr)
library(tidyr)
path <- "~/PJ/HJ/ZBP1_ZNA/01.counts/"
setwd(path)
  rm(list = ls())
  PS19 <- read.table("Zbp1_WT.feature_U_FC.txt", header = TRUE)
  head(PS19)
  PS19 <- unite(PS19, "geneid",
    c("Geneid", "Chr", "Start", "End"),
    sep = ",", remove = TRUE
  )
  head(PS19)
  Zbp1 <- read.table("Zbp1_PS19.feature_U_FC.txt", header = TRUE)
  head(Zbp1)
  Zbp1 <- unite(Zbp1, "geneid",
    c("Geneid", "Chr", "Start", "End"),
    sep = ",", remove = TRUE
  )
  head(Zbp1)
  exp <- cbind(
    PS19[
      ,
      c("geneid", "..00.bam.Zbp1_WTAligned.sortedByCoord.out.bam")
    ],
    Zbp1[
      ,
      c("..00.bam.Zbp1_PS19Aligned.sortedByCoord.out.bam")
    ]
  )
  head(exp)
  colnames(exp) <- c("Geneid", "WT", "Zbp1")
  exp <- exp[, c("Geneid", "WT", "Zbp1")]
  rownames(exp) <- exp[, 1]
  exp <- exp[, -1]
  head(exp)
  # 指定分组
  group <- rep(c("control", "treat"), each = 1)

  # 数据预处理
  # （1）构建 DGEList 对象
  dgelist <- DGEList(counts = exp, group = group)

  # （2）过滤 low count 数据，例如 CPM 标准化（推荐）
  keep <- rowSums(cpm(dgelist) > 1) >= 2
  dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]

  # （3）标准化
  dgelist_norm <- calcNormFactors(dgelist, method = "TMM")

  y_bcv <- dgelist
  # 因为本次的数据使用的是小鼠的数据，所以使用0.1
  bcv <- 0.01
  et <- exactTest(y_bcv, dispersion = bcv^2)
  results <- cbind(y_bcv$counts, et$table)

  # 将新生成的results数据框写成一个excel数据表
  write.csv(
    x = results,
    file = "Zbp1.edgeR.WTvsPS19.feature.csv", row.names = TRUE
  )
}
