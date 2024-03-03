rm(list = ls())
source("~/Desktop/Research/Rodu/MCCP/utils3.R")
load("~/Desktop/Research/Rodu/MCCP/Data/tar0.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/tar1.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/tar2.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/gp0.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/gp1.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/gp2.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/naive80.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/naive161.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/naive171.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/naive181.RData")
load("~/Desktop/Research/Rodu/MCCP/Data/naive191.RData")

## Extract the true change points for TAR and GP data
TAR <- GP <- vector("list", 3)
for(elem in 1:3) {
  TAR[[elem]] <- GP[[elem]]  <- vector("list", 100)
  for(sim in 1:100) {
    TAR[[elem]][[sim]] <- eval(parse(text = paste0("tar", elem - 1, "[[", sim, "]]$change_points")))
    GP[[elem]][[sim]] <- eval(parse(text = paste0("gp", elem - 1, "[[", sim, "]]$change_points")))
  }
}

## Extract the true change points for naive variance change data
ND <- vector("list", 5)
elems <- c(80, 161, 171, 181, 191)
for(elem in 1:5) {
  ND[[elem]] <- vector("list", 50)
  for(sim in 1:50) {
    ND[[elem]][[sim]] <- eval(parse(text = paste0("naive", elems[elem], "[[", sim, "]]$change_points")))
  }
}

## Visualize the distribution of change points for a given simulated data set
which <- TAR[[3]]
ggplot() +
  geom_histogram(aes(x = unlist(which)), binwidth = 40, center = 20, fill = "gray90", color = "black") +
  theme_bw()

## Load SCCP and MCCP results for TAR and GP simulations
for(process in c("t", "g")) {
  for(type in c("S", "M")) {
    main_folder <- ifelse(process == "t", "TAR", "GP")
    setwd(paste0("~/Desktop/Research/Rodu/MCCP/", main_folder, "/", type, "CCP/Results"))
    temp <- vector("list", 1200)
    for(sim in 1:1200) {
      if(file.exists(paste0(process, type, "CCP_", sim, ".RData"))) {
        load(paste0(process, type, "CCP_", sim, ".RData"))
        temp[[sim]] <- eval(parse(text = paste0(process, type, "CCP_", sim, "$estCPs")))
        rm(list = paste0(process, type, "CCP_", sim))
      } else {
        print(paste0("File does not exist: ", process, type, "CCP_", sim))
      }
    }
    assign(paste0(process, type, "CCP"), temp)
  }
}  

## Load MCCP results for naive variance change simulations
setwd("~/Desktop/Research/Rodu/MCCP/Naive/Results")
nMCCP <- vector("list", 1000)
for(sim in 1:1000) {
  if(file.exists(paste0("nMCCP_", sim, ".RData"))) {
    load(paste0("nMCCP_", sim, ".RData"))
    nMCCP[[sim]] <- eval(parse(text = paste0("nMCCP_", sim, "$estCPs")))
    rm(list = paste0("nMCCP_", sim))
  } else {
    print(paste0("File does not exist: nMCCP_", sim))
  }
}

## Load EDIV, KCP, NPMOJO, and SBS results
for(method in c("EDIV", "KCP", "NPMOJO", "SBS")) {
  setwd(paste0("~/Desktop/Research/Rodu/MCCP/", method))
  for(process in c("t", "g")) {
    for(ncp in 0:2) {
      if(file.exists(paste0(process, method, ncp, ".RData"))) {
        load(paste0(process, method, ncp, ".RData"))
      } else {
        print(paste0("File does not exist: ", process, method, ncp))
      }
    }
  }
}

# Clean and reconcile change point estimates that do not meet minimum spacing assumption
tEDIV0 <- cleanEDIV(tEDIV0)
tEDIV1 <- cleanEDIV(tEDIV1)
tEDIV2 <- cleanEDIV(tEDIV2)
gEDIV0 <- cleanEDIV(gEDIV0)
gEDIV1 <- cleanEDIV(gEDIV1)
gEDIV2 <- cleanEDIV(gEDIV2)
tKCP0 <- cleanKCP(tKCP0)
tKCP1 <- cleanKCP(tKCP1)
tKCP2 <- cleanKCP(tKCP2)
gKCP0 <- cleanKCP(gKCP0)
gKCP1 <- cleanKCP(gKCP1)
gKCP2 <- cleanKCP(gKCP2)
tNPMOJO0 <- cleanNPMOJO(tNPMOJO0)
tNPMOJO1 <- cleanNPMOJO(tNPMOJO1)
tNPMOJO2 <- cleanNPMOJO(tNPMOJO2)
gNPMOJO0 <- cleanNPMOJO(gNPMOJO0)
gNPMOJO1 <- cleanNPMOJO(gNPMOJO1)
gNPMOJO2 <- cleanNPMOJO(gNPMOJO2)
tSBS0 <- cleanSBS(tSBS0)
tSBS1 <- cleanSBS(tSBS1)
tSBS2 <- cleanSBS(tSBS2)
gSBS0 <- cleanSBS(gSBS0)
gSBS1 <- cleanSBS(gSBS1)
gSBS2 <- cleanSBS(gSBS2)

## Load ScanB and kCUSUM results
for(method in c("SB", "KC")) {
  setwd(paste0("~/Desktop/Research/Rodu/MCCP/", method))
  for(process in c("t", "g")) {
    for(ncp in 0:2) {
      for(ARL in 1:6) {
        if(file.exists(paste0(process, method, ncp, "_", ARL, ".csv"))) {
          assign("temp", split(read.csv(paste0(process, method, ncp, "_", ARL, ".csv"), header = FALSE), seq(nrow(read.csv(paste0(process, method, ncp, "_", ARL, ".csv"), header = FALSE)))))
          temp <- lapply(temp, function(x) x[(!is.na(x)) & (x != 0) & (x != 800)])
          assign(paste0(process, method, ncp, "_", ARL), temp)
        } else {
          print(paste0("File does not exist: ", process, method, ncp, ", ARL = ", ARL))
        }
      }
    }
  }
}

# Define radius for correct change point identification
radius <- 40

tON <- ONtable(TAR, nCP = "all", radius, "TAR")
round(tON, 3)
ggsave(plot = online_heatplot(tON), file = "~/Desktop/Research/Rodu/MCCP/Figures/TARonline.pdf", width = 10.5, height = 7.5)

tON0 <- ONtable(TAR, nCP = 0, radius, "TAR")
round(tON0, 3)
ggsave(plot = online_heatplot(tON0), file = "~/Desktop/Research/Rodu/MCCP/Figures/TARonline0.pdf", width = 10.5, height = 7.5)

tON1 <- ONtable(TAR, nCP = 1, radius, "TAR")
round(tON1, 3)
ggsave(plot = online_heatplot(tON1), file = "~/Desktop/Research/Rodu/MCCP/Figures/TARonline1.pdf", width = 10.5, height = 7.5)

tON2 <- ONtable(TAR, nCP = 2, radius, "TAR")
round(tON2, 3)
ggsave(plot = online_heatplot(tON2), file = "~/Desktop/Research/Rodu/MCCP/Figures/TARonline2.pdf", width = 10.5, height = 7.5)

gON <- ONtable(GP, nCP = "all", radius, "GP")
round(gON, 3)
ggsave(plot = online_heatplot(gON), file = "~/Desktop/Research/Rodu/MCCP/Figures/GPonline.pdf", width = 10.5, height = 7.5)

gON0 <- ONtable(GP, nCP = 0, radius, "GP")
round(gON0, 3)
ggsave(plot = online_heatplot(gON0), file = "~/Desktop/Research/Rodu/MCCP/Figures/GPonline0.pdf", width = 10.5, height = 7.5)

gON1 <- ONtable(GP, nCP = 1, radius, "GP")
round(gON1, 3)
ggsave(plot = online_heatplot(gON1), file = "~/Desktop/Research/Rodu/MCCP/Figures/GPonline1.pdf", width = 10.5, height = 7.5)

gON2 <- ONtable(GP, nCP = 2, radius, "GP")
round(gON2, 3)
ggsave(plot = online_heatplot(gON2), file = "~/Desktop/Research/Rodu/MCCP/Figures/GPonline2.pdf", width = 10.5, height = 7.5)

tOFF <- OFFtable(TAR, nCP = "all", radius, "TAR")
round(tOFF, 3)
ggsave(plot = offline_heatplot(tOFF), file = "~/Desktop/Research/Rodu/MCCP/Figures/TARoffline.pdf", width = 13, height = 7)

tOFF0 <- OFFtable(TAR, nCP = 0, radius, "TAR")
round(tOFF0, 3)
ggsave(plot = offline_heatplot(tOFF0), file = "~/Desktop/Research/Rodu/MCCP/Figures/TARoffline0.pdf", width = 13, height = 7)

tOFF1 <- OFFtable(TAR, nCP = 1, radius, "TAR")
round(tOFF1, 3)
ggsave(plot = offline_heatplot(tOFF1), file = "~/Desktop/Research/Rodu/MCCP/Figures/TARoffline1.pdf", width = 13, height = 7)

tOFF2 <- OFFtable(TAR, nCP = 2, radius, "TAR")
round(tOFF2, 3)
ggsave(plot = offline_heatplot(tOFF2), file = "~/Desktop/Research/Rodu/MCCP/Figures/TARoffline2.pdf", width = 13, height = 7)

gOFF <- OFFtable(GP, nCP = "all", radius, "GP")
round(gOFF, 3)
ggsave(plot = offline_heatplot(gOFF), file = "~/Desktop/Research/Rodu/MCCP/Figures/GPoffline.pdf", width = 13, height = 7)

gOFF0 <- OFFtable(GP, nCP = 0, radius, "GP")
round(gOFF0, 3)
ggsave(plot = offline_heatplot(gOFF0), file = "~/Desktop/Research/Rodu/MCCP/Figures/GPoffline0.pdf", width = 13, height = 7)

gOFF1 <- OFFtable(GP, nCP = 1, radius, "GP")
round(gOFF1, 3)
ggsave(plot = offline_heatplot(gOFF1), file = "~/Desktop/Research/Rodu/MCCP/Figures/GPoffline1.pdf", width = 13, height = 7)

gOFF2 <- OFFtable(GP, nCP = 2, radius, "GP")
round(gOFF2, 3)
ggsave(plot = offline_heatplot(gOFF2), file = "~/Desktop/Research/Rodu/MCCP/Figures/GPoffline2.pdf", width = 13, height = 7)


N <- Ntable(ND, radius)
round(N, 3)
quartz(file = "/Users/noahgade/Desktop/Research/Rodu/MCCP/Figures/Naive.pdf", type = "pdf", width = 12.8, height = 7)
grid.arrange(naive_heatplot(N))
dev.off()


## Load alternate combination results
setwd(paste0("~/Desktop/Research/Rodu/MCCP/NPMOJO"))
for(process in c("t", "g")) {
  for(ncp in 0:1) {
    if(file.exists(paste0(process, "CNPMOJO", ncp, ".RData"))) {
      load(paste0(process, "CNPMOJO", ncp, ".RData"))
    } else {
      print(paste0("File does not exist: ", process, "CNPMOJO", ncp))
    }
  }
}

tC <- Ctable(TAR, nCP = "all", radius, process = "TAR")
round(tC, 3)
ggsave(plot = Cheatplot(tC), file = "~/Desktop/Research/Rodu/MCCP/Figures/tC.pdf", width = 11, height = 3.8)

tC0 <- Ctable(TAR, nCP = 0, radius, process = "TAR")
round(tC0, 3)
ggsave(plot = Cheatplot(tC0), file = "~/Desktop/Research/Rodu/MCCP/Figures/tC0.pdf", width = 11, height = 3.8)

tC1 <- Ctable(TAR, nCP = 1, radius, process = "TAR")
round(tC1, 3)
ggsave(plot = Cheatplot(tC1), file = "~/Desktop/Research/Rodu/MCCP/Figures/tC1.pdf", width = 11, height = 3.8)

gC <- Ctable(GP, nCP = "all", radius, process = "GP")
round(gC, 3)
ggsave(plot = Cheatplot(gC), file = "~/Desktop/Research/Rodu/MCCP/Figures/gC.pdf", width = 11, height = 3.8)

gC0 <- Ctable(GP, nCP = 0, radius, process = "GP")
round(gC0, 3)
ggsave(plot = Cheatplot(gC0), file = "~/Desktop/Research/Rodu/MCCP/Figures/gC0.pdf", width = 11, height = 3.8)

gC1 <- Ctable(GP, nCP = 1, radius, process = "GP")
round(gC1, 3)
ggsave(plot = Cheatplot(gC1), file = "~/Desktop/Research/Rodu/MCCP/Figures/gC1.pdf", width = 11, height = 3.8)
