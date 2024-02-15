options(max.print = 1000)
options(stringsAsFactors = FALSE)
options(scipen = 999)

## description
# small random forest classifier trained on demo dataset (medulloblastoma subgroups, hovestadt 2013)
# 107 samples belonging to 4 tumor classes, 48 differentially methylated probes

library("randomForest")


## load and validate data
betas <- read.csv("/Users/volker/DropboxPartners/Benfatto_Figures/github/401_2013_1126_MOESM2_ESM.csv", row.names = 1)  # demo dataset
str(betas)
dim(betas)

y <- factor(sapply(strsplit(rownames(betas), "_"), head, 1), c("WNT", "SHH", "Group3", "Group4"))
table(y)

heatmap(t(as.matrix(betas)), scale="n", col=colorRampPalette(c("blue", "white", "red"))(100), useRaster=TRUE)  # as expected


## train classifer
set.seed(123)
rf.pred <- randomForest(betas, y, ntree = 1000, do.trace = TRUE, sampsize = rep(10, 4), keep.inbag = TRUE)   # downsample each class to 10 samples, inbag needed for later!
str(rf.pred)

save(betas, file="betas.RData")
save(rf.pred, file="rf.RData")
