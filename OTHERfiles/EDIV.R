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
  a01 <- e.divisive(data$data[161:640,], sig.lvl = 0.01)$estimates
  a01 <- a01[(a01 != 1) & (a01 != 481)] + 160 - 1
  b01 <- e.divisive(data$data[161:640,], sig.lvl = 0.01, min.size = 160)$estimates
  b01 <- b01[(b01 != 1) & (b01 != 481)] + 160 - 1
  a05 <- e.divisive(data$data[161:640,], sig.lvl = 0.05)$estimates
  a05 <- a05[(a05 != 1) & (a05 != 481)] + 160 - 1
  b05 <- e.divisive(data$data[161:640,], sig.lvl = 0.05, min.size = 160)$estimates
  b05 <- b05[(b05 != 1) & (b05 != 481)] + 160 - 1
  a10 <- e.divisive(data$data[161:640,], sig.lvl = 0.10)$estimates
  a10 <- a10[(a10 != 1) & (a10 != 481)] + 160 - 1
  b10 <- e.divisive(data$data[161:640,], sig.lvl = 0.10, min.size = 160)$estimates
  b10 <- b10[(b10 != 1) & (b10 != 481)] + 160 - 1
  return(list(a01 = a01, b01 = b01, a05 = a05, b05 = b05, a10 = a10, b10 = b10))
}

load("~/Desktop/Research/Rodu/MCCP/Data/tar0.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/tar1.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/tar2.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/gp0.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/gp1.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/gp2.RData")

data <- tar0
tEDIV0 <- vector("list", length(data))
pb <- progress_bar$new(total = length(data))
for(i in 1:length(data)) {
  tEDIV0[[i]] <- method(data[[i]])
  pb$tick()
}
save(tEDIV0, file = "~/Desktop/Research/Rodu/MCCP/EDIV/tEDIV0.RData")

data <- tar1
tEDIV1 <- vector("list", length(data))
pb <- progress_bar$new(total = length(data))
for(i in 1:length(data)) {
  tEDIV1[[i]] <- method(data[[i]])
  pb$tick()
}
save(tEDIV1, file = "~/Desktop/Research/Rodu/MCCP/EDIV/tEDIV1.RData")

data <- tar2
tEDIV2 <- vector("list", length(data))
pb <- progress_bar$new(total = length(data))
for(i in 1:length(data)) {
  tEDIV2[[i]] <- method(data[[i]])
  pb$tick()
}
save(tEDIV2, file = "~/Desktop/Research/Rodu/MCCP/EDIV/tEDIV2.RData")

data <- gp0
gEDIV0 <- vector("list", length(data))
pb <- progress_bar$new(total = length(data))
for(i in 1:length(data)) {
  gEDIV0[[i]] <- method(data[[i]])
  pb$tick()
}
save(gEDIV0, file = "~/Desktop/Research/Rodu/MCCP/EDIV/gEDIV0.RData")

data <- gp1
gEDIV1 <- vector("list", length(data))
pb <- progress_bar$new(total = length(data))
for(i in 1:length(data)) {
  gEDIV1[[i]] <- method(data[[i]])
  pb$tick()
}
save(gEDIV1, file = "~/Desktop/Research/Rodu/MCCP/EDIV/gEDIV1.RData")

data <- gp2
gEDIV2 <- vector("list", length(data))
pb <- progress_bar$new(total = length(data))
for(i in 1:length(data)) {
  gEDIV2[[i]] <- method(data[[i]])
  pb$tick()
}
save(gEDIV2, file = "~/Desktop/Research/Rodu/MCCP/EDIV/gEDIV2.RData")
