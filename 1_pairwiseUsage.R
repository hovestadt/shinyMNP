options(max.print = 1000)
options(stringsAsFactors = FALSE)
options(scipen = 999)

#library("parallel")
library("randomForest")
library("data.table")
library("iterpc")

## script to extract summary statistics from all trees of a random forest
# pairwise depth between all samples: at which point in trees are samples split
# pairwise jaccard between all samples: how much of the path in trees is shared
# pairwise probe usage between all classes: how often is a probe selected to split


## load classifier
f <- commandArgs(trailingOnly = TRUE)[1]
# f <- "varselRFlistv11b2.mtry655.seed14541013.RData"  # outer/probe selection classifier (all probes)
f <- "rfpred.v11b2.prox.RData"  # inner classifier (10,000 probes)

# which forests from list (outer classifier)
# if(length(commandArgs(trailingOnly = TRUE)) == 2) {
#   rr <- commandArgs(trailingOnly = TRUE)[2]
# } else {
#   rr <- "1-43"``
# }
# rr <- as.numeric(strsplit(rr, "-")[[1]])
rr <- 1

# # quit if single forest and already existing
# Sys.sleep(runif(1, 0, 10))
# if(length(rr) == 1 & file.exists(paste0(f, ".pairComb-", rr,".RData"))) {
#   q()
# } else {
#   save(f, file=paste0(f, ".pairComb-", rr, ".RData"))  # dummy
# }


## load predictions, if exist, or generate
if(file.exists(paste0(f, ".predNode.RData"))) {
  # includes rf and predictions
  load(paste0(f, ".predNode.RData"))
} else {
  load(f)
  
  # probe selection classifier:
  # 43* ~10k probes, mtry: 655 (sqrt of all used probes)
  # 100 trees each
  # final classifier:
  # 10k probes, mtry: 100 (sqrt of 10k)
  # 10000 trees
  
  # put into list if not already (rfpred.v11b2.prox.RData)
  if(!is.element("rflist", ls())) {
    rflist <- list(rf.pred)
    rm("rf.pred")
  }
  
  load("MNPv11b2_2801samples.20160325.betas.ba.RData")
  gc()
  rf.pred.node.list <- lapply(rflist, function(rf.pred) {
    attr(predict(rf.pred, betas[, rownames(rf.pred$importance)], nodes = TRUE), "nodes")
  })
  rm("betas")
  gc()
  
  save(list=c("rflist", "rf.pred.node.list"), file=paste0(f, ".predNode.RData"))
}

str(rflist)
str(rf.pred.node.list)
head(rf.pred.node.list[[1]][, 1])  # terminal node in tree 1


## functions
extractPaths <- function(k.tree) {
  # returns a list of vectors for each terminal node in the supplied tree
  # values correspond to row numbers in tree
  k.tree.term <- which(k.tree$status == -1)
  k.tree.paths <- lapply(k.tree.term, function(path) {
    #path <- 267  # terminal node
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
  # second and third matrix contain the probe and wheather it is hypo/hypermethylated
  k.tree.paths.combn <- getall(iterpc(length(k.tree.paths), 2, replace=TRUE, ordered=FALSE))  # only do for upper diag, the while loop is somewhat slow for all comparisons. fill lower diag later
  #dim(k.tree.paths.combn)
  
  k.tree.paths.combn.split <- do.call(rbind,
                                      mapply(k.tree.paths[k.tree.paths.combn[, 1]],
                                             k.tree.paths[k.tree.paths.combn[, 2]],
                                             FUN = function(path1, path2) {
                                               #i <- 1217
                                               #path1 <- k.tree.paths[[k.tree.paths.combn[i, 1]]]
                                               #path2 <- k.tree.paths[[k.tree.paths.combn[i, 2]]]
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
  
  # all pairwise sample combinations (2801x2801 = 7845601)
  sp <- getall(iterpc(length(s.node), 2, replace=TRUE, ordered=TRUE))  # full matrix, to keep things simpler. same for every tree, move outside of loop?
  sp.inbag <-  s.inbag[sp[, 1]] * s.inbag[sp[, 2]]
  
  # fill data.table with in-bag sample combinations (8x8 x 91x91 = 529984)
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


## loop over defined forests in list (all if not defined)
for(r in seq(head(rr, 1), tail(rr, 1))) {
  r <- 1  # forest
  message(paste0("* forest ", r))
  
  rf.pred <- rflist[[r]]
  rf.pred.node <- rf.pred.node.list[[r]]
  
  ## init objects
  sp_inb.depth <- matrix(0, ncol=length(rf.pred$y), nrow=length(rf.pred$y))  # SamplePairwise, two dimensional
  dim(sp_inb.depth)
  sp_inb.jaccard <- sp_inb.depth
  sp_inb.n <- sp_inb.depth
  
  cp_inb.probe <- array(0, dim = c(length(rf.pred$classes), length(rf.pred$classes), length(rf.pred$forest$ncat)))  # ClassPairwise, three dimensional (incl probe)
  dim(cp_inb.probe)
  # there is no n; probe usage is summed up, using -/+1 for hypo/hypermethylated 
  
  ## loop over trees
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
    
    sp_inb.depth[k.tree.paths.split.pair$sample[, 1:2]] <- sp_inb.depth[k.tree.paths.split.pair$sample[, 1:2]] + k.tree.paths.split.pair$sample[, 3]
    sp_inb.jaccard[k.tree.paths.split.pair$sample[, 1:2]] <- sp_inb.jaccard[k.tree.paths.split.pair$sample[, 1:2]] + k.tree.paths.split.pair$sample[, 4]
    sp_inb.n[k.tree.paths.split.pair$sample[, 1:2]] <- sp_inb.n[k.tree.paths.split.pair$sample[, 1:2]] + k.tree.paths.split.pair$sample[, 5]
    
    cp_inb.probe[k.tree.paths.split.pair$class[, 1:3]] <- cp_inb.probe[k.tree.paths.split.pair$class[, 1:3]] + k.tree.paths.split.pair$class[, 4]
      
    #invisible(gc())
  }
  
  rflist.pair <- list(depth=sp_inb.depth, jaccard=sp_inb.jaccard, n=sp_inb.n, probe=cp_inb.probe,
                      sample_names=rownames(rf.pred.node), class_names=rf.pred$classes, probe_names=names(rf.pred$forest$ncat))
  
  save(rflist.pair, file=paste0(f, ".pairComb-", r,".RData"))
  invisible(gc())
}

str(rflist.pair)

