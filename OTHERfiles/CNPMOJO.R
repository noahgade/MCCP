rm(list = ls())
source("~/Desktop/Research/Rodu/MCCP/utils3.R")
sourceCpp("~/Desktop/Research/Rodu/MCCP/MCCP.cpp")
source("~/Desktop/Research/Rodu/MCCP/utils2.R")

method <- function(data) {
  a01 <- cnp.mojo.ESM(data$data, washL = 40, trainL = 120, threshold = 0.01)
  a05 <- cnp.mojo.ESM(data$data, washL = 40, trainL = 120, threshold = 0.05)
  a10 <- cnp.mojo.ESM(data$data, washL = 40, trainL = 120, threshold = 0.10)
  return(list(a01 = a01, a05 = a05, a10 = a10))
}

load("~/Desktop/Research/Rodu/MCCP/Data/tar0.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/tar1.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/gp0.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/gp1.RData")

data <- tar0
tCNPMOJO0 <- vector("list", length(data))
pb <- progress_bar$new(total = length(data))
for(i in 1:length(data)) {
  tCNPMOJO0[[i]] <- method(data[[i]])
  pb$tick()
}
save(tCNPMOJO0, file = "~/Desktop/Research/Rodu/MCCP/NPMOJO/tCNPMOJO0.RData")

data <- tar1
tCNPMOJO1 <- vector("list", length(data))
pb <- progress_bar$new(total = length(data))
for(i in 1:length(data)) {
  tCNPMOJO1[[i]] <- method(data[[i]])
  pb$tick()
}
save(tCNPMOJO1, file = "~/Desktop/Research/Rodu/MCCP/NPMOJO/tCNPMOJO1.RData")

data <- gp0
gCNPMOJO0 <- vector("list", length(data))
pb <- progress_bar$new(total = length(data))
for(i in 1:length(data)) {
  gCNPMOJO0[[i]] <- method(data[[i]])
  pb$tick()
}
save(gCNPMOJO0, file = "~/Desktop/Research/Rodu/MCCP/NPMOJO/gCNPMOJO0.RData")

data <- gp1
gCNPMOJO1 <- vector("list", length(data))
pb <- progress_bar$new(total = length(data))
for(i in 1:length(data)) {
  gCNPMOJO1[[i]] <- method(data[[i]])
  pb$tick()
}
save(gCNPMOJO1, file = "~/Desktop/Research/Rodu/MCCP/NPMOJO/gCNPMOJO1.RData")
