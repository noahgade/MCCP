rm(list = ls())
library(parallel)
numCores <- detectCores()
options(mc.cores = numCores)
library(progress)
library(tidyverse)
library(matrixStats)
library(ecp)
library(hdbinseg)
library(CptNonPar)

method <- function(data) {
  r101 <- sbs.alg(t(data$data[161:640,]), cp.type = 1, q = 0.01)$ecp + 160
  r105 <- sbs.alg(t(data$data[161:640,]), cp.type = 1, q = 0.05)$ecp + 160
  r110 <- sbs.alg(t(data$data[161:640,]), cp.type = 1, q = 0.10)$ecp + 160
  r201 <- sbs.alg(t(data$data[161:640,]), cp.type = 2, q = 0.01)$ecp + 160
  r205 <- sbs.alg(t(data$data[161:640,]), cp.type = 2, q = 0.05)$ecp + 160
  r210 <- sbs.alg(t(data$data[161:640,]), cp.type = 2, q = 0.10)$ecp + 160
  return(list(r101 = r101, r105 = r105, r110 = r110,
              r201 = r201, r205 = r205, r210 = r210))
}

load("~/Desktop/Research/Rodu/MCCP/Data/tar0.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/tar1.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/tar2.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/gp0.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/gp1.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/gp2.RData")

data <- tar0
tSBS0 <- vector("list", length(data))
pb <- progress_bar$new(total = length(data))
for(i in 1:length(data)) {
  tSBS0[[i]] <- method(data[[i]])
  pb$tick()
}
save(tSBS0, file = "~/Desktop/Research/Rodu/MCCP/SBS/tSBS0.RData")

data <- tar1
tSBS1 <- vector("list", length(data))
pb <- progress_bar$new(total = length(data))
for(i in 1:length(data)) {
  tSBS1[[i]] <- method(data[[i]])
  pb$tick()
}
save(tSBS1, file = "~/Desktop/Research/Rodu/MCCP/SBS/tSBS1.RData")

data <- tar2
tSBS2 <- vector("list", length(data))
pb <- progress_bar$new(total = length(data))
for(i in 1:length(data)) {
  tSBS2[[i]] <- method(data[[i]])
  pb$tick()
}
save(tSBS2, file = "~/Desktop/Research/Rodu/MCCP/SBS/tSBS2.RData")

data <- gp0
gSBS0 <- vector("list", length(data))
pb <- progress_bar$new(total = length(data))
for(i in 1:length(data)) {
  gSBS0[[i]] <- method(data[[i]])
  pb$tick()
}
save(gSBS0, file = "~/Desktop/Research/Rodu/MCCP/SBS/gSBS0.RData")

data <- gp1
gSBS1 <- vector("list", length(data))
pb <- progress_bar$new(total = length(data))
for(i in 1:length(data)) {
  gSBS1[[i]] <- method(data[[i]])
  pb$tick()
}
save(gSBS1, file = "~/Desktop/Research/Rodu/MCCP/SBS/gSBS1.RData")

data <- gp2
gSBS2 <- vector("list", length(data))
pb <- progress_bar$new(total = length(data))
for(i in 1:length(data)) {
  gSBS2[[i]] <- method(data[[i]])
  pb$tick()
}
save(gSBS2, file = "~/Desktop/Research/Rodu/MCCP/SBS/gSBS2.RData")

