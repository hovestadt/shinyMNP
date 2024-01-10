options(max.print = 1000)
options(stringsAsFactors = FALSE)
options(scipen = 999)

## description
# script to extract the following summary statistics from all trees of a random forest:
# pairwise depth between all samples: at which point in trees are samples split (not used in publication)
# pairwise jaccard between all samples: how much of the path in trees is shared (also not used)
# pairwise probe usage between all classes: how often is a probe selected to split (used throughout)
#
# volker hovestadt, salvatore benfatto
# dkfz, dfci
# 2018-2024

library("randomForest")
library("data.table")
library("iterpc")

source("shinyMNP_functions.R")


## load rf classifier
f <- "rf.v11b2.RData"  # brain tumor classifier (10,000 probes), replace with other rf object of interest
load(f, verbose = TRUE)
rf.pred
# Type of random forest: classification
# Number of trees: 10000
# No. of variables tried at each split: 100


## define terminal node for every reference sample (in-bag and out-of-bag samples)
load("betas.v11b2.RData", verbose = TRUE)  # exact reference that was used to train the rf classifier
dim(betas)  # sample x probe
# [1]   2801 428799

rf.pred.node <- attr(predict(rf.pred, betas[, rownames(rf.pred$importance)], nodes = TRUE), "nodes")
dim(rf.pred.node)  # sample x tree
# [1]  2801 10000
head(rf.pred.node[, 1])  # terminal node for every sample in tree 1

rm(betas)  # not needed anymore
gc()


## define functions
extractPaths <- function(k.tree) {
  # returns a list of vectors for each terminal node in the supplied tree
  # values correspond to row numbers in tree
  k.tree.term <- which(k.tree$status == -1)
  k.tree.paths <- lapply(k.tree.term, function(path) {
    # path <- 267  # terminal node
    while(tail(path, 1) > 1) {
      w1 <- k.tree[ "left daughter"] == tail(path, 1)
      if(any(w1)) {
        path <- c(path, which(w1))
      } else {
        w2 <- k.tree[, "right daughter"] == tail(path, 1)
        path <- c(path, which(w2))
      }
    }
    rev(path)  # starts at 1 (root), ends at terminal node
  })
  k.tree.paths
}

comparePaths <- function(k.tree.paths, k.tree) {
  # comparison of pairwise combinations of paths
  # first matrix in list gives the split depth of paths, diagonal contains depth of terminal nodes
  # second and third matrix contain the probe and whether it is hypo/hypermethylated
  k.tree.paths.combn <- getall(iterpc(length(k.tree.paths), 2, replace=TRUE, ordered=FALSE))  # only do for upper diag, the while loop is somewhat slow for all comparisons. fill lower diag later
  # dim(k.tree.paths.combn)
  
  k.tree.paths.combn.split <- do.call(rbind,
                                      mapply(k.tree.paths[k.tree.paths.combn[, 1]],
                                             k.tree.paths[k.tree.paths.combn[, 2]],
                                             FUN = function(path1, path2) {
                                               # i <- 1217
                                               # path1 <- k.tree.paths[[k.tree.paths.combn[i, 1]]]
                                               # path2 <- k.tree.paths[[k.tree.paths.combn[i, 2]]]
                                               d <- 0  # split depth
                                               try(while(path1[d+1] == path2[d+1]) d <- d+1, silent = TRUE)  # catch error for exact same paths: terminal node depth
                                               d.probe <- path1[d]  # split probe (same in both paths), or terminal node
                                               d.hyper1 <- ifelse(is.element(path1[d+1], k.tree[, "right daughter"]), 1, -1)  # -1 if path1 is hypomehtylated at split node, +1 of hyper, check for terminal nodes below
                                               c(d-1, d.probe, d.hyper1)
                                             }, SIMPLIFY = FALSE, USE.NAMES = FALSE))
  k.tree.paths.combn.split[is.na(k.tree[k.tree.paths.combn.split[, 2], "split var"]), 3] <- 0  # set hyper to 0 for terminal nodes
  
  k.tree.paths.depth <- k.tree.paths.probe <- k.tree.paths.hyper <- matrix(NA, ncol=length(k.tree.paths), nrow=length(k.tree.paths))
  k.tree.paths.depth[k.tree.paths.combn] <- k.tree.paths.combn.split[, 1]
  k.tree.paths.depth[lower.tri(k.tree.paths.depth)] <- t(k.tree.paths.depth)[lower.tri(k.tree.paths.depth)]  # also fill lower.tri, pairwise samples could be either way
  k.tree.paths.probe[k.tree.paths.combn] <- match(as.vector(k.tree[, "split var"]), names(rf.pred$forest$ncat))[k.tree.paths.combn.split[, 2]]  # use probe no from full classifier
  k.tree.paths.probe[lower.tri(k.tree.paths.probe)] <- t(k.tree.paths.probe)[lower.tri(k.tree.paths.probe)]
  k.tree.paths.hyper[k.tree.paths.combn] <- k.tree.paths.combn.split[, 3]
  k.tree.paths.hyper[lower.tri(k.tree.paths.hyper)] <- -t(k.tree.paths.hyper)[lower.tri(k.tree.paths.hyper)]  # change sign!
  
  list(depth=k.tree.paths.depth, probe=k.tree.paths.probe, hyper=k.tree.paths.hyper)
}

depth2jaccard <- function(k.tree.paths.split.depth) {
  # calculates jaccard index of shared decisions over all decisions from split depth matrix
  matrix(nrow = nrow(k.tree.paths.split.depth), as.vector(k.tree.paths.split.depth) /
    (rep(diag(k.tree.paths.split.depth), each = nrow(k.tree.paths.split.depth)) +
       rep(diag(k.tree.paths.split.depth), nrow(k.tree.paths.split.depth)) -
       as.vector(k.tree.paths.split.depth)))
}

pairComb <- function(k.tree.paths.split, k.tree.paths.split.jaccard, k.tree, s.node, s.inbag) {
  # extrapolate pairwise paths/terminal node combinations to all pairwise in-bag sample combinations
  # if a sample is in-bag multiple times, this is taken into account
  # usage of probes is aggregated to class combinations (output too large for sample combinations)
  s.node <- unname(s.node)
  s.inbag <- unname(s.inbag)
  s.term <- match(s.node, which(k.tree$status == -1))
  
  # all pairwise sample combinations (2,801x2,801 = 7,845,601)
  sp <- getall(iterpc(length(s.node), 2, replace=TRUE, ordered=TRUE))  # full matrix, to keep things simpler. same for every tree, move outside of loop?
  sp.inbag <-  s.inbag[sp[, 1]] * s.inbag[sp[, 2]]
  
  # fill data.table with in-bag sample combinations (8x8 x 91x91 = 529,984)
  sp.rep <- data.table(sampleA = rep(sp[, 1], sp.inbag), sampleB = rep(sp[, 2], sp.inbag))  # rep by inbag, need to use aggregate:sum later!
  sp.rep$classA <- as.numeric(rf.pred$y)[sp.rep$sampleA]
  sp.rep$classB <- as.numeric(rf.pred$y)[sp.rep$sampleB]
  sp.rep$n <- 1
  sp.rep.term <- cbind(s.term[sp.rep$sampleA], s.term[sp.rep$sampleB])
  sp.rep$depth <- k.tree.paths.split$depth[sp.rep.term]
  sp.rep$jaccard <- k.tree.paths.split.jaccard[sp.rep.term]
  sp.rep$probe <- k.tree.paths.split$probe[sp.rep.term]
  sp.rep$hyper <- k.tree.paths.split$hyper[sp.rep.term]  # sum should always be zero
  
  # aggregate by sample and class+probe
  sp.rep.sample <- sp.rep[, lapply(.SD, sum), by=c("sampleA", "sampleB"), .SDcols=c("depth", "jaccard", "n")] 
  sp.rep.class <- sp.rep[, lapply(.SD, sum), by=c("classA", "classB", "probe"), .SDcols=c("hyper")]  # samples ending up in the same terminal node (including the same sample compared to itself) have value 0, no need to filter
  
  return(list(sample = as.matrix(sp.rep.sample), class = as.matrix(sp.rep.class[!is.na(sp.rep.class$probe), ])))  # remove terminal nodes, convert to matrix
}


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


## save output file
rflist.pair <- list(depth=sp_inb.depth, jaccard=sp_inb.jaccard, n=sp_inb.n, probe=cp_inb.probe,
                    sample_names=rownames(rf.pred.node), class_names=rf.pred$classes, probe_names=names(rf.pred$forest$ncat))
str(rflist.pair)

save(rflist.pair, file=paste0(f, ".pairComb.RData"))
gc()


## some simple summaries
# probe usage between classes
dim(cp_inb.probe)  # class1 x class2 x probe
# [1]    91    91 10000

# total probe usage
sum(cp_inb.probe)
range(cp_inb.probe)
sum(abs(cp_inb.probe))

# most used probes
range(apply(abs(cp_inb.probe), 3, sum))
plot(sort(apply(abs(cp_inb.probe), 3, sum)), type="l", log="y")

which.max(apply(abs(cp_inb.probe), 3, sum))
image(cp_inb.probe[, , 48], col=colorRampPalette(c("blue", "white", "red"))(100))  # most highly used probe
image(cp_inb.probe[, , 1234], col=colorRampPalette(c("blue", "white", "red"))(100))   # random probe




