## description
# functions used to analyze trees of random forest classifier
#
# volker hovestadt, salvatore benfatto
# dkfz, dfci
# 2018-2024


## define functions
message("extractPaths()")
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

message("comparePaths()")
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

message("depth2jaccard()")
depth2jaccard <- function(k.tree.paths.split.depth) {
  # calculates jaccard index of shared decisions over all decisions from split depth matrix
  matrix(nrow = nrow(k.tree.paths.split.depth), as.vector(k.tree.paths.split.depth) /
           (rep(diag(k.tree.paths.split.depth), each = nrow(k.tree.paths.split.depth)) +
              rep(diag(k.tree.paths.split.depth), nrow(k.tree.paths.split.depth)) -
              as.vector(k.tree.paths.split.depth)))
}

message("pairComb()")
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
