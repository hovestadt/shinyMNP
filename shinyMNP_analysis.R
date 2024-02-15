options(max.print = 1000)
options(stringsAsFactors = FALSE)
options(scipen = 999)

## description
# script to extract the following summary statistics from all trees of a random forest:
# pairwise depth between all samples: at which point in trees are samples split (not used in publication)
# pairwise jaccard between all samples: how much of the path in trees is shared (also not used)
# pairwise probe usage between all classes: how often is a probe selected to split (used throughout publication)
#
# volker hovestadt, salvatore benfatto
# dkfz, dfci
# 2018-2024

library("randomForest")
library("data.table")
library("iterpc")

source("shinyMNP_functions.R")


## load rf classifier
f <- "rf.RData"  # demo classifier (medulloblastoma subgroups, hovestadt 2013), replace with other rf object of interest
load(f, verbose = TRUE)
rf.pred
# Type of random forest: classification
# Number of trees: 1000
# No. of variables tried at each split: 6


## define terminal node for every reference sample (in-bag and out-of-bag samples)
load("betas.RData", verbose = TRUE)  # exact reference that was used to train the rf classifier
dim(betas)  # sample x probe
# [1] 107  48

rf.pred.node <- attr(predict(rf.pred, betas[, rownames(rf.pred$importance)], nodes = TRUE), "nodes")
dim(rf.pred.node)  # sample x tree
# [1]  107 1000
head(rf.pred.node[, 1])  # terminal node for every sample in tree 1

rm(betas)  # not needed anymore
gc()


## init objects
sp_inb.depth <- matrix(0, ncol=length(rf.pred$y), nrow=length(rf.pred$y))  # SamplePairwise, two dimensional
dim(sp_inb.depth)  # sample x sample
sp_inb.jaccard <- sp_inb.depth
sp_inb.n <- sp_inb.depth

cp_inb.probe <- array(0, dim = c(length(rf.pred$classes), length(rf.pred$classes), length(rf.pred$forest$ncat)))  # ClassPairwise, three dimensional (incl probe)
dim(cp_inb.probe)  # class x class
# there is no n; probe usage is summed up, using -/+1 for hypo/hypermethylated 


## loop over trees in forest (~1sec/tree)
for(k in seq(ncol(rf.pred.node))) {
  # k <- 1  # tree
  
  ## get tree
  message(paste0("+ tree ", k))
  k.tree <- getTree(rf.pred, k, labelVar=TRUE)
  # str(k.tree)
  
  ## calculating path for every terminal node in tree, quite fast
  message("- extracting paths")
  k.tree.paths <- extractPaths(k.tree)
  # str(k.tree.paths)
  
  ## compare all pairs of paths
  message("- comparing paths")
  k.tree.paths.split <- comparePaths(k.tree.paths, k.tree)
  # str(k.tree.paths.split)
  
  ## calculating jaccard index
  message("- calculating jaccard")
  k.tree.paths.split.jaccard <- depth2jaccard(k.tree.paths.split$depth)
  # str(k.tree.paths.split.jaccard)
  
  ## compare all pairwise in-bag sample combinations
  message("- comparing all sample combinations")
  k.tree.paths.split.pair <- pairComb(k.tree.paths.split, k.tree.paths.split.jaccard, k.tree, rf.pred.node[, k], rf.pred$inbag[, k])
  #str(k.tree.paths.split.pair)
  
  ## add values to preexisting matrices
  message("- adding values")
  sp_inb.depth[k.tree.paths.split.pair$sample[, 1:2]] <- sp_inb.depth[k.tree.paths.split.pair$sample[, 1:2]] + k.tree.paths.split.pair$sample[, 3]  # needs to be divided by n later
  sp_inb.jaccard[k.tree.paths.split.pair$sample[, 1:2]] <- sp_inb.jaccard[k.tree.paths.split.pair$sample[, 1:2]] + k.tree.paths.split.pair$sample[, 4]
  sp_inb.n[k.tree.paths.split.pair$sample[, 1:2]] <- sp_inb.n[k.tree.paths.split.pair$sample[, 1:2]] + k.tree.paths.split.pair$sample[, 5]
  
  cp_inb.probe[k.tree.paths.split.pair$class[, 1:3]] <- cp_inb.probe[k.tree.paths.split.pair$class[, 1:3]] + k.tree.paths.split.pair$class[, 4]
  invisible(gc())
}


## save output files
rflist.pair <- list(depth=sp_inb.depth, jaccard=sp_inb.jaccard, n=sp_inb.n, probe=cp_inb.probe,
                    sample_names=rownames(rf.pred.node), class_names=rf.pred$classes, probe_names=names(rf.pred$forest$ncat))
str(rflist.pair)

save(rflist.pair, file=paste0(f, ".pairComb.RData"))


cp_inb.probe.df <- as.data.frame(t(apply(rflist.pair$probe, 3, unlist)))
rownames(cp_inb.probe.df) <- rflist.pair$probe_names
colnames(cp_inb.probe.df) <- paste0(rep(rflist.pair$class_names, 4), " vs ", rep(rflist.pair$class_names, each=4))

write.table(cp_inb.probe.df, file=paste0(f, ".pairComb.txt"), sep="\t", quote=FALSE)


## some simple summaries
# probe usage between classes
dim(cp_inb.probe)  # class1 x class2 x probe
# [1]  4  4 48

# total probe usage
sum(cp_inb.probe)   # has to be zero
range(cp_inb.probe)
sum(abs(cp_inb.probe))

# most used probes overall
range(apply(abs(cp_inb.probe), 3, sum))
plot(sort(apply(abs(cp_inb.probe), 3, sum)), type="l", log="y")

which.max(apply(abs(cp_inb.probe), 3, sum))
image(cp_inb.probe[, , 40], col=colorRampPalette(c("blue", "white", "red"))(100))  # most highly used probe

# probes distinguishing classes
plot(cp_inb.probe[3, 4,])  # group 3 from group 4
plot(colSums(cp_inb.probe[3, ,]))  # group 3 from all others

