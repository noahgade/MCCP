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
  a01 <- np.mojo.multilag(data$data, G = 160, lags = 0:2, alpha = 0.01)
  a05 <- np.mojo.multilag(data$data, G = 160, lags = 0:2, alpha = 0.05)
  a10 <- np.mojo.multilag(data$data, G = 160, lags = 0:2, alpha = 0.10)
  return(list(a01 = a01, a05 = a05, a10 = a10))
}

load("~/Desktop/Research/Rodu/MCCP/Data/tar0.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/tar1.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/tar2.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/gp0.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/gp1.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/gp2.RData")

data <- tar0
tNPMOJO0 <- vector("list", length(data))
pb <- progress_bar$new(total = length(data))
for(i in 1:length(data)) {
  tNPMOJO0[[i]] <- method(data[[i]])
  pb$tick()
}
save(tNPMOJO0, file = "~/Desktop/Research/Rodu/MCCP/NPMOJO/tNPMOJO0.RData")

data <- tar1
tNPMOJO1 <- vector("list", length(data))
pb <- progress_bar$new(total = length(data))
for(i in 1:length(data)) {
  tNPMOJO1[[i]] <- method(data[[i]])
  pb$tick()
}
save(tNPMOJO1, file = "~/Desktop/Research/Rodu/MCCP/NPMOJO/tNPMOJO1.RData")

data <- tar2
tNPMOJO2 <- vector("list", length(data))
pb <- progress_bar$new(total = length(data))
for(i in 1:length(data)) {
  tNPMOJO2[[i]] <- method(data[[i]])
  pb$tick()
}
save(tNPMOJO2, file = "~/Desktop/Research/Rodu/MCCP/NPMOJO/tNPMOJO2.RData")

data <- gp0
gNPMOJO0 <- vector("list", length(data))
pb <- progress_bar$new(total = length(data))
for(i in 1:length(data)) {
  gNPMOJO0[[i]] <- method(data[[i]])
  pb$tick()
}
save(gNPMOJO0, file = "~/Desktop/Research/Rodu/MCCP/NPMOJO/gNPMOJO0.RData")

data <- gp1
gNPMOJO1 <- vector("list", length(data))
pb <- progress_bar$new(total = length(data))
for(i in 1:length(data)) {
  gNPMOJO1[[i]] <- method(data[[i]])
  pb$tick()
}
save(gNPMOJO1, file = "~/Desktop/Research/Rodu/MCCP/NPMOJO/gNPMOJO1.RData")

data <- gp2
gNPMOJO2 <- vector("list", length(data))
pb <- progress_bar$new(total = length(data))
for(i in 1:length(data)) {
  gNPMOJO2[[i]] <- method(data[[i]])
  pb$tick()
}
save(gNPMOJO2, file = "~/Desktop/Research/Rodu/MCCP/NPMOJO/gNPMOJO2.RData")
