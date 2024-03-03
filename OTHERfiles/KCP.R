rm(list = ls())
library(parallel)
numCores <- detectCores()
options(mc.cores = numCores)
library(tidyverse)
library(matrixStats)
library(ecp)
library(hdbinseg)

method <- function(data) {
  r02 <- sapply(1:5, function(x) kcpa(as.matrix(data$data[161:640,]), L = x, C = 2) - 1 + 160, simplify = F)
  r20 <- sapply(1:5, function(x) kcpa(as.matrix(data$data[161:640,]), L = x, C = 20) - 1 + 160, simplify = F)
  r40 <- sapply(1:5, function(x) kcpa(as.matrix(data$data[161:640,]), L = x, C = 40) - 1 + 160, simplify = F)
  r60 <- sapply(1:5, function(x) kcpa(as.matrix(data$data[161:640,]), L = x, C = 60) - 1 + 160, simplify = F)
  r80 <- sapply(1:5, function(x) kcpa(as.matrix(data$data[161:640,]), L = x, C = 80) - 1 + 160, simplify = F)
  r100 <- sapply(1:5, function(x) kcpa(as.matrix(data$data[161:640,]), L = x, C = 100) - 1 + 160, simplify = F)
  return(list(r02 = r02, r20 = r20, r40 = r40, r60 = r60, r80 = r80, r100 = r100))
}

load("tar0.RData")
load("tar1.RData")
load("tar2.RData")
load("gp0.RData")
load("gp1.RData")
load("gp2.RData")

tKCP0 <- vector("list", length(tar0))
for(i in 1:length(tar0)) {
  tKCP0[[i]] <- method(tar0[[i]])
}
save(tKCP0, file = "~/home/ndg5e/MCCP/tKCP0.RData")

tKCP1 <- vector("list", length(tar1))
for(i in 1:length(tar1)) {
  tKCP1[[i]] <- method(tar1[[i]])
}
save(tKCP1, file = "~/home/ndg5e/MCCP/tKCP1.RData")

tKCP2 <- vector("list", length(tar2))
for(i in 1:length(tar2)) {
  tKCP2[[i]] <- method(tar2[[i]])
}
save(tKCP2, file = "~/home/ndg5e/MCCP/tKCP2.RData")

gKCP0 <- vector("list", length(gp0))
for(i in 1:length(gp0)) {
  gKCP0[[i]] <- method(gp0[[i]])
}
save(gKCP0, file = "~/home/ndg5e/MCCP/gKCP0.RData")

gKCP1 <- vector("list", length(gp1))
for(i in 1:length(gp1)) {
  gKCP1[[i]] <- method(gp1[[i]])
}
save(gKCP1, file = "~/home/ndg5e/MCCP/gKCP1.RData")

gKCP2 <- vector("list", length(gp2))
for(i in 1:length(gp2)) {
  gKCP2[[i]] <- method(gp2[[i]])
}
save(gKCP2, file = "~/home/ndg5e/MCCP/gKCP2.RData")
