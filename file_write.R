rm(list = ls())
library(data.table)
library(stringr)
n_sims <- 100
n_simsN <- 50
CPs <- c(0, 1, 2)
error <- c(TRUE, FALSE)
threshold <- c(0.05, 0.10)
naive <- c(80, 161, 171, 181, 191)

set.seed(224)
tSCCPseeds <- sample(0:10000, length(CPs) * length(threshold) * length(error) * n_sims, replace = TRUE)
tMCCPseeds <- sample(0:10000, length(CPs) * length(threshold) * length(error) * n_sims, replace = TRUE)
gSCCPseeds <- sample(0:10000, length(CPs) * length(threshold) * length(error) * n_sims, replace = TRUE)
gMCCPseeds <- sample(0:10000, length(CPs) * length(threshold) * length(error) * n_sims, replace = TRUE)
nMCCPseeds <- sample(0:10000, length(naive) * length(threshold) * length(error) * n_simsN, replace = TRUE)

for(change in 1:length(CPs)) {
  for(control in 1:length(error)) {
    for(alpha in 1:length(threshold)) {
      for(sim in 1:n_sims) {
        count <- (change - 1) * length(error) * length(threshold) * n_sims + (control - 1) * length(threshold) * n_sims + (alpha - 1) * n_sims + sim
        filename <- paste0("~/Desktop/Research/Rodu/MCCP/TAR/SCCP/tSCCP_", count, ".R")
        fileConn <- file(filename)
        writeLines(c("numCores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK')) - 1",
                     "options(mc.cores = numCores)", 
                     "library(parallel)",
                     "library(tidyverse)",
                     "library(matrixStats)",
                     "library(abind)",
                     "library(Rcpp)",
                     "library(RcppArmadillo)",
                     "sourceCpp('/home/ndg5e/MCCP/MCCP.cpp')",
                     paste0("load('/home/ndg5e/MCCP/tar", CPs[change], ".RData')"),
                     paste0("dat <- tar", CPs[change], "[[", sim, "]]$data"),
                     paste0("set.seed(", tSCCPseeds[count], ")"),
                     paste0("tSCCP_", count, " <- SCCP(dat, washL = 40, trainL = 120, threshold = ", threshold[alpha], ", FWER = ", error[control], ")"),
                     paste0("save(tSCCP_", count, ", file = '/home/ndg5e/MCCP/tSCCPresults/tSCCP_", count, ".RData')"),
                     "print('COMPLETE')"),
                   fileConn) 
        close(fileConn)
      }
    }
  }
}

for(change in 1:length(CPs)) {
  for(control in 1:length(error)) {
    for(alpha in 1:length(threshold)) {
      for(sim in 1:n_sims) {
        count <- (change - 1) * length(error) * length(threshold) * n_sims + (control - 1) * length(threshold) * n_sims + (alpha - 1) * n_sims + sim
        filename <- paste0("~/Desktop/Research/Rodu/MCCP/TAR/MCCP/tMCCP_", count, ".R")
        fileConn <- file(filename)
        writeLines(c("numCores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK')) - 1",
                     "options(mc.cores = numCores)", 
                     "library(parallel)",
                     "library(tidyverse)",
                     "library(matrixStats)",
                     "library(abind)",
                     "library(Rcpp)",
                     "library(RcppArmadillo)",
                     "sourceCpp('/home/ndg5e/MCCP/MCCP.cpp')",
                     paste0("load('/home/ndg5e/MCCP/tar", CPs[change], ".RData')"),
                     paste0("dat <- tar", CPs[change], "[[", sim, "]]$data"),
                     paste0("set.seed(", tMCCPseeds[count], ")"),
                     paste0("tMCCP_", count, " <- MCCP(dat, washL = 40, trainL = 120, threshold = ", threshold[alpha], ", FWER = ", error[control], ")"),
                     paste0("save(tMCCP_", count, ", file = '/home/ndg5e/MCCP/tMCCPresults/tMCCP_", count, ".RData')"),
                     "print('COMPLETE')"),
                   fileConn) 
        close(fileConn)
      }
    }
  }
}

for(change in 1:length(CPs)) {
  for(control in 1:length(error)) {
    for(alpha in 1:length(threshold)) {
      for(sim in 1:n_sims) {
        count <- (change - 1) * length(error) * length(threshold) * n_sims + (control - 1) * length(threshold) * n_sims + (alpha - 1) * n_sims + sim
        filename <- paste0("~/Desktop/Research/Rodu/MCCP/GP/SCCP/gSCCP_", count, ".R")
        fileConn <- file(filename)
        writeLines(c("numCores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK')) - 1",
                     "options(mc.cores = numCores)", 
                     "library(parallel)",
                     "library(tidyverse)",
                     "library(matrixStats)",
                     "library(abind)",
                     "library(Rcpp)",
                     "library(RcppArmadillo)",
                     "sourceCpp('/home/ndg5e/MCCP/MCCP.cpp')",
                     paste0("load('/home/ndg5e/MCCP/gp", CPs[change], ".RData')"),
                     paste0("dat <- gp", CPs[change], "[[", sim, "]]$data"),
                     paste0("set.seed(", gSCCPseeds[count], ")"),
                     paste0("gSCCP_", count, " <- SCCP(dat, washL = 40, trainL = 120, threshold = ", threshold[alpha], ", FWER = ", error[control], ")"),
                     paste0("save(gSCCP_", count, ", file = '/home/ndg5e/MCCP/gSCCPresults/gSCCP_", count, ".RData')"),
                     "print('COMPLETE')"),
                   fileConn) 
        close(fileConn)
      }
    }
  }
}

for(change in 1:length(CPs)) {
  for(control in 1:length(error)) {
    for(alpha in 1:length(threshold)) {
      for(sim in 1:n_sims) {
        count <- (change - 1) * length(error) * length(threshold) * n_sims + (control - 1) * length(threshold) * n_sims + (alpha - 1) * n_sims + sim
        filename <- paste0("~/Desktop/Research/Rodu/MCCP/GP/MCCP/gMCCP_", count, ".R")
        fileConn <- file(filename)
        writeLines(c("numCores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK')) - 1",
                     "options(mc.cores = numCores)", 
                     "library(parallel)",
                     "library(tidyverse)",
                     "library(matrixStats)",
                     "library(abind)",
                     "library(Rcpp)",
                     "library(RcppArmadillo)",
                     "sourceCpp('/home/ndg5e/MCCP/MCCP.cpp')",
                     paste0("load('/home/ndg5e/MCCP/gp", CPs[change], ".RData')"),
                     paste0("dat <- gp", CPs[change], "[[", sim, "]]$data"),
                     paste0("set.seed(", gMCCPseeds[count], ")"),
                     paste0("gMCCP_", count, " <- MCCP(dat, washL = 40, trainL = 120, threshold = ", threshold[alpha], ", FWER = ", error[control], ")"),
                     paste0("save(gMCCP_", count, ", file = '/home/ndg5e/MCCP/gMCCPresults/gMCCP_", count, ".RData')"),
                     "print('COMPLETE')"),
                   fileConn) 
        close(fileConn)
      }
    }
  }
}

for(spacing in 1:length(naive)) {
  for(control in 1:length(error)) {
    for(alpha in 1:length(threshold)) {
      for(sim in 1:n_simsN) {
        count <- (spacing - 1) * length(error) * length(threshold) * n_simsN + (control - 1) * length(threshold) * n_simsN + (alpha - 1) * n_simsN + sim
        filename <- paste0("~/Desktop/Research/Rodu/MCCP/Naive/nMCCP_", count, ".R")
        fileConn <- file(filename)
        writeLines(c("numCores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK')) - 1",
                     "options(mc.cores = numCores)",
                     "library(parallel)",
                     "library(tidyverse)",
                     "library(matrixStats)",
                     "library(abind)",
                     "library(Rcpp)",
                     "library(RcppArmadillo)",
                     "sourceCpp('/home/ndg5e/MCCP/MCCP.cpp')",
                     paste0("load('/home/ndg5e/MCCP/naive", naive[spacing], ".RData')"),
                     paste0("dat <- naive", naive[spacing], "[[", sim, "]]$data"),
                     paste0("set.seed(", nMCCPseeds[count], ")"),
                     paste0("nMCCP_", count, " <- MCCP(dat, washL = 40, trainL = 120, threshold = ", threshold[alpha], ", FWER = ", error[control], ")"),
                     paste0("save(nMCCP_", count, ", file = '/home/ndg5e/MCCP/nMCCPresults/nMCCP_", count, ".RData')"),
                     "print('COMPLETE')"),
                   fileConn)
        close(fileConn)
      }
    }
  }
}
