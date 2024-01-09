options(max.print = 1000)
options(stringsAsFactors = FALSE)
options(scipen = 999)

load("../MNPv11b2_2801samples.20160325.betas.ba.RData")
rm("v11_anno")
gc()

betas_names <- sort(colnames(betas))
length(betas_names)
rm("betas")
gc()

sp_inb.depth.sum <- matrix(0, ncol=2801, nrow=2801)
sp_inb.jaccard.sum <- sp_inb.depth.sum
sp_inb.n.sum <- sp_inb.depth.sum

cp_inb.probe.sum <- array(0, dim = c(91, 91, length(betas_names)))
gc()

lf <- list.files(pattern = "pairComb")
length(lf)

for(f in lf) {
  # f <- lf[1]
  message(f)
  load(f)
  
  sp_inb.depth.sum <- sp_inb.depth.sum + rflist.pair$depth
  sp_inb.jaccard.sum <- sp_inb.jaccard.sum + rflist.pair$jaccard
  sp_inb.n.sum <- sp_inb.n.sum + rflist.pair$n
  
  i <- match(rflist.pair$probe_names, betas_names)
  cp_inb.probe.sum[, , i] <- cp_inb.probe.sum[, , i] + rflist.pair$probe
  
  invisible(gc())
}

save(list=c("sp_inb.depth.sum", "sp_inb.jaccard.sum", "sp_inb.n.sum"), file="all_pairComb_jaccard.vh20160727.RData")
save(list=c("cp_inb.probe.sum", "betas_names"), file="all_pairComb_probe.vh20160727.RData")
