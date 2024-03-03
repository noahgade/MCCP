rm(list = ls())
library(R.matlab)
load("~/Desktop/Research/Rodu/MCCP/Data/tar0.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/tar1.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/tar2.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/gp0.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/gp1.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/gp2.RData")

convert <- function(data_list) {
  out <- array(dim = c(ncol(data_list[[1]]$data), nrow(data_list[[1]]$data), length(data_list)))
  for(i in 1:length(data_list)) {
    out[,,i] <- t(data_list[[i]]$data)
  }
  return(out)
}

writeMat("~/Desktop/Research/Rodu/MCCP/Data/tar0.mat", tar0 = convert(tar0))
writeMat("~/Desktop/Research/Rodu/MCCP/Data/tar1.mat", tar1 = convert(tar1))
writeMat("~/Desktop/Research/Rodu/MCCP/Data/tar2.mat", tar2 = convert(tar2))
writeMat("~/Desktop/Research/Rodu/MCCP/Data/gp0.mat", gp0 = convert(gp0))
writeMat("~/Desktop/Research/Rodu/MCCP/Data/gp1.mat", gp1 = convert(gp1))
writeMat("~/Desktop/Research/Rodu/MCCP/Data/gp2.mat", gp2 = convert(gp2))