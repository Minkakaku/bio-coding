# install.packages("stringr")
require(stringr)
d <- "/home/molab/hongfan/project/YXX/results/2"
setwd(d)
df <- read.csv("diff_paper2.csv",
    stringsAsFactors = F,
    header = T
)
df$ensembl_id <- unlist(str_split(df$ID, "[.]", simplify = T))[, 1]

write.csv(df, file = "noedition.csv")
