rm(list = ls())
source("~/Desktop/Research/Rodu/MCCP/utils.R")

n_sims <- 100
n_simsN <- 50
washL <- 40
trainL <- 120

# sample(0:10000, 1)
# [1] 7415
set.seed(7415)
TARMat <- vector("list", 3 * n_sims)
pb_TARMat <- progress_bar$new(total = length(TARMat))
for(sim in 1:(3 * n_sims)) {
  TARMat[[sim]] <- genTARmat(0.8)
  pb_TARMat$tick()
}
save(TARMat, file = "~/Desktop/Research/Rodu/MCCP/Data/TARMat.RData")

# sample(0:10000, 1)
# [1] 5833
set.seed(5833)
tar0 <- vector("list", n_sims)
tar1 <- vector("list", n_sims)
tar2 <- vector("list", n_sims)
for(sim in 1:n_sims) {
  tar0[[sim]] <- genTAR(length = 800, n_changes = 0, min_spacing = washL + trainL + 1, TARmat = TARMat[[sim]])
  tar1[[sim]] <- genTAR(length = 800, n_changes = 1, min_spacing = washL + trainL + 1, TARmat = TARMat[[n_sims + sim]])
  tar2[[sim]] <- genTAR(length = 800, n_changes = 2, min_spacing = washL + trainL + 1, TARmat = TARMat[[2 * n_sims + sim]])
}
save(tar0, file = "~/Desktop/Research/Rodu/MCCP/Data/tar0.RData")
save(tar1, file = "~/Desktop/Research/Rodu/MCCP/Data/tar1.RData")
save(tar2, file = "~/Desktop/Research/Rodu/MCCP/Data/tar2.RData")

# sample(0:10000, 1)
# [1] 8340
set.seed(8340)
gp0 <- vector("list", n_sims)
gp1 <- vector("list", n_sims)
gp2 <- vector("list", n_sims)
for(sim in 1:n_sims) {
  gp0[[sim]] <- genGP(length = 800, n_changes = 0, min_spacing = washL + trainL + 1)
  gp1[[sim]] <- genGP(length = 800, n_changes = 1, min_spacing = washL + trainL + 1)
  gp2[[sim]] <- genGP(length = 800, n_changes = 2, min_spacing = washL + trainL + 1)
  print(sim)
}
save(gp0, file = "~/Desktop/Research/Rodu/MCCP/Data/gp0.RData")
save(gp1, file = "~/Desktop/Research/Rodu/MCCP/Data/gp1.RData")
save(gp2, file = "~/Desktop/Research/Rodu/MCCP/Data/gp2.RData")

# sample(0:10000, 1)
# [1] 6529
set.seed(6529)
naive80 <- vector("list", n_simsN)
naive161 <- vector("list", n_simsN)
naive171 <- vector("list", n_simsN)
naive181 <- vector("list", n_simsN)
naive191 <- vector("list", n_simsN)
for(sim in 1:n_simsN) {
  naive80[[sim]] <- genNaive(length = 800, spacing = 80)
  naive161[[sim]] <- genNaive(length = 800, spacing = 161)
  naive171[[sim]] <- genNaive(length = 800, spacing = 171)
  naive181[[sim]] <- genNaive(length = 800, spacing = 181)
  naive191[[sim]] <- genNaive(length = 800, spacing = 191)
}
save(naive80, file = "~/Desktop/Research/Rodu/MCCP/Data/naive80.RData")
save(naive161, file = "~/Desktop/Research/Rodu/MCCP/Data/naive161.RData")
save(naive171, file = "~/Desktop/Research/Rodu/MCCP/Data/naive171.RData")
save(naive181, file = "~/Desktop/Research/Rodu/MCCP/Data/naive181.RData")
save(naive191, file = "~/Desktop/Research/Rodu/MCCP/Data/naive191.RData")
