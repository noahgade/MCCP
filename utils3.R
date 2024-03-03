source("~/Desktop/Research/Rodu/MCCP/utils.R")
sourceCpp("~/Desktop/Research/Rodu/MCCP/MCCP.cpp")

score <- function(ListEST, ListTRUE, radius = 1, npts = 480) {
  if(class(ListEST) == "numeric") {
    ListEST <- list(ListEST)
    ListTRUE <- list(ListTRUE)
    print("Warning: change point vectors converted to lists")
  }
  len <- min(length(ListEST), length(ListTRUE))
  if(length(ListEST) != length(ListTRUE)) {
    ListEST <- ListEST[1:len]
    ListTRUE <- ListTRUE[1:len]
    print("Warning: unequal lengths - minimum list length chosen for analysis")
  }
  TPR <- FPR <- FDR <- FOR <- MCC <- rep(0, len)
  for(elem in 1:len) {
    lenEST <- length(ListEST[[elem]])
    lenTRUE <- length(ListTRUE[[elem]])
    if((lenEST == 0) & (lenTRUE == 0)) {
      TPR[elem] <- 1    # Scenario without true positives
      FPR[elem] <- 0
      FDR[elem] <- 0
      FOR[elem] <- 0
    } else if((lenEST == 0) & (lenTRUE > 0)) {
      TPR[elem] <- 0
      FPR[elem] <- 0
      FDR[elem] <- 0
      FOR[elem] <- lenTRUE / npts
    } else if((lenEST > 0) & (lenTRUE == 0)) {
      TPR[elem] <- 1    # Scenario without true positives
      FPR[elem] <- lenEST / npts
      FDR[elem] <- 1
      FOR[elem] <- 0
    } else {
      nTP <- 0
      for(cp in 1:lenTRUE) {
        dists <- abs(ListEST[[elem]] - ListTRUE[[elem]][cp])
        if(min(dists) <= radius) {
          nTP <- nTP + 1
        }
      }
      nTP <- min(nTP, lenEST)
      TPR[elem] <- nTP / lenTRUE
      FPR[elem] <- (lenEST - nTP) / (npts - lenTRUE)
      FDR[elem] <- (lenEST - nTP) / lenEST
      FOR[elem] <- (lenTRUE - nTP) / (npts - lenEST)
    }
    MCC[elem] <- sqrt(TPR[elem] * (1 - FPR[elem]) * (1 - FDR[elem]) * (1 - FOR[elem]))  - sqrt((1 - TPR[elem]) * FPR[elem] * FDR[elem] * FOR[elem])
  }
  return(c(TPR = mean(TPR), FPR = mean(FPR), FDR = mean(FDR), FOR = mean(FOR), MCC = mean(MCC)))
}

rreconcile <- function(FWD, BWD, radius, L) {
  out <- Reconcile(matrix(FWD), matrix(BWD), radius, L)
  out <- out[(out > 0) & (out < L)]
  return(out)
}

cleanEDIV <- function(EDIV) {
  len <- length(EDIV)
  for(i in 1:len) {
    EDIV[[i]]$a01 <- rreconcile(EDIV[[i]]$a01, EDIV[[i]]$a01, 161, 800)
    EDIV[[i]]$a05 <- rreconcile(EDIV[[i]]$a05, EDIV[[i]]$a05, 161, 800)
    EDIV[[i]]$a10 <- rreconcile(EDIV[[i]]$a10, EDIV[[i]]$a10, 161, 800)
  }
  return(EDIV)
}

cleanKCP <- function(KCP) {
  len1 <- length(KCP)
  len2 <- length(KCP[[1]])
  len3 <- length(KCP[[1]][[1]])
  for(i in 1:len1) {
    for(j in 1:len2) {
      for(k in 1:len3) {
        KCP[[i]][[j]][[k]] <- KCP[[i]][[j]][[k]][(KCP[[i]][[j]][[k]] > 160) & (KCP[[i]][[j]][[k]] < 640)]
      }
    }
  }
  output <- vector("list", length = len1)
  for(i in 1:len1) {
    output[[i]] <- vector("list", length = len2)
    for(j in 1:len2) {
      r <- vector("list", length = len3)
      for(k in 1:len3) {
        r[[k]] <- rreconcile(KCP[[i]][[j]][[k]], KCP[[i]][[j]][[k]], 161, 800)
      }
      output[[i]][[j]] <- r[[min(which(lapply(r, function(x) length(x)) == 2), 5)]]
    }
  }
  return(output)
}

cleanNPMOJO <- function(NPMOJO) {
  len <- length(NPMOJO)
  out <- vector("list", length = len)
  for(i in 1:len) {
    out[[i]] <- vector("list", length = 3)
    nr1 <- nrow(NPMOJO[[i]]$a01$cpts)
    nr2 <- nrow(NPMOJO[[i]]$a05$cpts)
    nr3 <- nrow(NPMOJO[[i]]$a10$cpts)
    out[[i]][[1]] <- as.vector(NPMOJO[[i]]$a01$cpts[,1])
    out[[i]][[2]] <- as.vector(NPMOJO[[i]]$a05$cpts[,1])
    out[[i]][[3]] <- as.vector(NPMOJO[[i]]$a10$cpts[,1])
  }
  return(out)
}

cleanSBS <- function(SBS) {
  len <- length(SBS)
  for(i in 1:len) {
    SBS[[i]]$r101 <- rreconcile(SBS[[i]]$r101, SBS[[i]]$r101, 161, 800)
    SBS[[i]]$r105 <- rreconcile(SBS[[i]]$r105, SBS[[i]]$r105, 161, 800)
    SBS[[i]]$r110 <- rreconcile(SBS[[i]]$r110, SBS[[i]]$r110, 161, 800)
    SBS[[i]]$r201 <- rreconcile(SBS[[i]]$r201, SBS[[i]]$r201, 161, 800)
    SBS[[i]]$r205 <- rreconcile(SBS[[i]]$r205, SBS[[i]]$r205, 161, 800)
    SBS[[i]]$r210 <- rreconcile(SBS[[i]]$r210, SBS[[i]]$r210, 161, 800)
  }
  return(SBS)
}

ONtable <- function(trueCP, nCP = "all", radius, process = "TAR") {
  class <- ifelse(process == "TAR", "t", "g")
  cp <- ifelse(nCP == "all", "", nCP)
  SCCP <- eval(parse(text = paste0(class, "SCCP")))
  for(method in c("SB", "KC")) {
    for(ARL in 1:6) {
      assign(paste0(method, ARL), eval(parse(text = paste0("c(", class, method, "0_", ARL, ", ", class, method, "1_", ARL, ", ", class, method, "2_", ARL, ")"))))
    }
  }
  SCCPA <- c(SCCP[1:100], SCCP[401:500], SCCP[801:900])
  SCCPB <- c(SCCP[101:200], SCCP[501:600], SCCP[901:1000])
  SCCPC <- c(SCCP[201:300], SCCP[601:700], SCCP[1001:1100])
  SCCPD <- c(SCCP[301:400], SCCP[701:800], SCCP[1101:1200])
  tCP <- c(trueCP[[1]], trueCP[[2]], trueCP[[3]])
  if(nCP == 0) {
    range <- 1:100
  } else if(nCP == 1) {
    range <- 101:200
  } else if(nCP == 2) {
    range <- 201:300
  } else {
    range <- 1:300
  }
  tON <- as.table(rbind(score(SCCPA[range], tCP[range], radius, npts = 480),
                        score(SCCPB[range], tCP[range], radius, npts = 480),
                        score(SCCPC[range], tCP[range], radius, npts = 480),
                        score(SCCPD[range], tCP[range], radius, npts = 480),
                        score(SB1[range], tCP[range], radius, npts = 480),
                        score(SB2[range], tCP[range], radius, npts = 480),
                        score(SB3[range], tCP[range], radius, npts = 480),
                        score(SB4[range], tCP[range], radius, npts = 480),
                        score(SB5[range], tCP[range], radius, npts = 480),
                        score(SB6[range], tCP[range], radius, npts = 480),
                        score(KC1[range], tCP[range], radius, npts = 480),
                        score(KC2[range], tCP[range], radius, npts = 480),
                        score(KC3[range], tCP[range], radius, npts = 480),
                        score(KC4[range], tCP[range], radius, npts = 480),
                        score(KC5[range], tCP[range], radius, npts = 480),
                        score(KC6[range], tCP[range], radius, npts = 480),
                        c(1, 0, 0, 0, 1)))
  dimnames(tON) <- list(Method = c("SCCP[FWER]_q = 0.05", "SCCP[FWER]_q = 0.10",
                                  "SCCP[FDR]_q = 0.05", "SCCP[FDR]_q = 0.10",
                                  "ScanB_ARL = 1e+01", "ScanB_ARL = 1e+02", "ScanB_ARL = 1e+03", 
                                  "ScanB_ARL = 1e+04", "ScanB_ARL = 1e+05", "ScanB_ARL = 1e+06",
                                  "kCUSUM_ARL = 1e+01", "kCUSUM_ARL = 1e+02", "kCUSUM_ARL = 1e+03",
                                  "kCUSUM_ARL = 1e+04", "kCUSUM_ARL = 1e+05", "kCUSUM_ARL = 1e+06",
                                  paste0("Oracle_", process, cp)),
                        Metric = c("TPR", "FPR", "FDR", "FOR", "MCC"))
return(tON)
}

OFFtable <- function(trueCP, nCP = "all", radius, process = "TAR") {
  class <- ifelse(process == "TAR", "t", "g")
  cp <- ifelse(nCP == "all", "", nCP)
  MCCP <- eval(parse(text = paste0(class, "MCCP")))
  for(method in c("EDIV", "KCP", "NPMOJO", "SBS")) {
    assign(paste0(method), eval(parse(text = paste0("c(", class, method, "0", ", ", class, method, "1", ", ", class, method, "2", ")"))))
  }
  MCCPA <- c(MCCP[1:100], MCCP[401:500], MCCP[801:900])
  MCCPB <- c(MCCP[101:200], MCCP[501:600], MCCP[901:1000])
  MCCPC <- c(MCCP[201:300], MCCP[601:700], MCCP[1001:1100])
  MCCPD <- c(MCCP[301:400], MCCP[701:800], MCCP[1101:1200])
  tCP <- c(trueCP[[1]], trueCP[[2]], trueCP[[3]])
  if(nCP == 0) {
    range <- 1:100
  } else if(nCP == 1) {
    range <- 101:200
  } else if(nCP == 2) {
    range <- 201:300
  } else {
    range <- 1:300
  }
  tOF <- as.table(rbind(score(MCCPA[range], tCP[range], radius, npts = 480),
                        score(MCCPB[range], tCP[range], radius, npts = 480),
                        score(MCCPC[range], tCP[range], radius, npts = 480),
                        score(MCCPD[range], tCP[range], radius, npts = 480),
                        # score(lapply(EDIV[range], function(x) x$a01), tCP[range], radius, npts = 480),
                        # score(lapply(EDIV[range], function(x) x$a05), tCP[range], radius, npts = 480),
                        # score(lapply(EDIV[range], function(x) x$a10), tCP[range], radius, npts = 480),
                        score(lapply(EDIV[range], function(x) x$b01), tCP[range], radius, npts = 480),
                        score(lapply(EDIV[range], function(x) x$b05), tCP[range], radius, npts = 480),
                        score(lapply(EDIV[range], function(x) x$b10), tCP[range], radius, npts = 480),
                        score(lapply(KCP[range], function(x) x[[6]]), tCP[range], radius, npts = 480),
                        score(lapply(KCP[range], function(x) x[[5]]), tCP[range], radius, npts = 480),
                        score(lapply(KCP[range], function(x) x[[4]]), tCP[range], radius, npts = 480),
                        score(lapply(KCP[range], function(x) x[[3]]), tCP[range], radius, npts = 480),
                        score(lapply(KCP[range], function(x) x[[2]]), tCP[range], radius, npts = 480),
                        score(lapply(KCP[range], function(x) x[[1]]), tCP[range], radius, npts = 480),
                        score(lapply(NPMOJO[range], function(x) x[[1]]), tCP[range], radius, npts = 480),
                        score(lapply(NPMOJO[range], function(x) x[[2]]), tCP[range], radius, npts = 480),
                        score(lapply(NPMOJO[range], function(x) x[[3]]), tCP[range], radius, npts = 480),
                        score(lapply(SBS[range], function(x) x$r101), tCP[range], radius, npts = 480),
                        score(lapply(SBS[range], function(x) x$r105), tCP[range], radius, npts = 480),
                        score(lapply(SBS[range], function(x) x$r110), tCP[range], radius, npts = 480),
                        score(lapply(SBS[range], function(x) x$r201), tCP[range], radius, npts = 480),
                        score(lapply(SBS[range], function(x) x$r205), tCP[range], radius, npts = 480),
                        score(lapply(SBS[range], function(x) x$r210), tCP[range], radius, npts = 480),
                        c(1, 0, 0, 0, 1)))
  dimnames(tOF) <- list(Method = c("MCCP[FWER]_q = 0.05", "MCCP[FWER]_q = 0.10",
                                   "MCCP[FDR]_q = 0.05", "MCCP[FDR]_q = 0.10",
                                   # "EDiv[MS = 30]_q = 0.01", "EDiv[MS = 30]_q = 0.05", "EDiv[MS = 30]_q = 0.10",
                                   "EDiv[MS = 160]_q = 0.01", "EDiv[MS = 160]_q = 0.05", "EDiv[MS = 160]_q = 0.10",
                                   "KCP_C = 100", "KCP_C = 80", 
                                   "KCP_C = 60", "KCP_C = 40", 
                                   "KCP_C = 20", "KCP_C = 2",
                                   "NP-MOJO[G = 160]_q = 0.01", "NP-MOJO[G = 160]_q = 0.05", "NP-MOJO[G = 160]_q = 0.10",
                                   "SBS1_q = 0.01", "SBS1_q = 0.05", "SBS1_q = 0.10",
                                   "SBS2_q = 0.01", "SBS2_q = 0.05", "SBS2_q = 0.10",
                                   paste0("Oracle_", process, cp)),
                        Metric = c("TPR", "FPR", "FDR", "FOR", "MCC"))
  return(tOF)
}

Ntable <- function(trueCP, radius) {
  nMCCPA <- c(nMCCP[1:50], nMCCP[201:250], nMCCP[401:450], nMCCP[601:650], nMCCP[801:850])
  nMCCPB <- c(nMCCP[51:100], nMCCP[251:300], nMCCP[451:500], nMCCP[651:700], nMCCP[851:900])
  nMCCPC <- c(nMCCP[101:150], nMCCP[301:350], nMCCP[501:550], nMCCP[701:750], nMCCP[901:950])
  nMCCPD <- c(nMCCP[151:200], nMCCP[351:400], nMCCP[551:600], nMCCP[751:800], nMCCP[951:1000])
  tCP <- c(trueCP[[1]], trueCP[[2]], trueCP[[3]], trueCP[[4]], trueCP[[5]])

  tOF <- as.table(rbind(score(nMCCP[1:50], trueCP[[1]], radius, npts = 480),
                        score(nMCCP[201:250], trueCP[[2]], radius, npts = 480),
                        score(nMCCP[401:450], trueCP[[3]], radius, npts = 480),
                        score(nMCCP[601:650], trueCP[[4]], radius, npts = 480),
                        score(nMCCP[801:850], trueCP[[5]], radius, npts = 480),
                        score(nMCCP[51:100], trueCP[[1]], radius, npts = 480),
                        score(nMCCP[251:300], trueCP[[2]], radius, npts = 480),
                        score(nMCCP[451:500], trueCP[[3]], radius, npts = 480),
                        score(nMCCP[651:700], trueCP[[4]], radius, npts = 480),
                        score(nMCCP[851:900], trueCP[[5]], radius, npts = 480),
                        score(nMCCP[101:150], trueCP[[1]], radius, npts = 480),
                        score(nMCCP[301:350], trueCP[[2]], radius, npts = 480),
                        score(nMCCP[501:550], trueCP[[3]], radius, npts = 480),
                        score(nMCCP[701:750], trueCP[[4]], radius, npts = 480),
                        score(nMCCP[901:950], trueCP[[5]], radius, npts = 480),
                        score(nMCCP[151:200], trueCP[[1]], radius, npts = 480),
                        score(nMCCP[351:400], trueCP[[2]], radius, npts = 480),
                        score(nMCCP[551:600], trueCP[[3]], radius, npts = 480),
                        score(nMCCP[751:800], trueCP[[4]], radius, npts = 480),
                        score(nMCCP[951:1000], trueCP[[5]], radius, npts = 480),
                        c(1, 0, 0, 0, 1)))
  dimnames(tOF) <- list(Method = c("FWER[q = 0.05]_\u03b3 = 80", "FWER[q = 0.05]_\u03b3 = 161", "FWER[q = 0.05]_\u03b3 = 171", "FWER[q = 0.05]_\u03b3 = 181", "FWER[q = 0.05]_\u03b3 = 191",
                                   "FWER[q = 0.10]_\u03b3 = 80", "FWER[q = 0.10]_\u03b3 = 161", "FWER[q = 0.10]_\u03b3 = 171", "FWER[q = 0.10]_\u03b3 = 181", "FWER[q = 0.10]_\u03b3 = 191",
                                   "FDR[q = 0.05]_\u03b3 = 80", "FDR[q = 0.05]_\u03b3 = 161", "FDR[q = 0.05]_\u03b3 = 171", "FDR[q = 0.05]_\u03b3 = 181", "FDR[q = 0.05]_\u03b3 = 191",
                                   "FDR[q = 0.10]_\u03b3 = 80", "FDR[q = 0.10]_\u03b3 = 161", "FDR[q = 0.10]_\u03b3 = 171", "FDR[q = 0.10]_\u03b3 = 181", "FDR[q = 0.10]_\u03b3 = 191",
                                   "Oracle_Naive"),
                        Metric = c("TPR", "FPR", "FDR", "FOR", "MCC"))
  return(tOF)
}

Ctable <- function(trueCP, nCP = "all", radius, process = "TAR") {
  class <- ifelse(process == "TAR", "t", "g")
  cp <- ifelse(nCP == "all", "", nCP)
  MCCP <- eval(parse(text = paste0(class, "MCCP")))
  CNPMOJO <- eval(parse(text = paste0("c(", class, "CNPMOJO", 0, ", ", class, "CNPMOJO", 1, ")")))
  NPMOJO <- eval(parse(text = paste0("c(", class, "NPMOJO", 0, ", ", class, "NPMOJO", 1, ")")))

  MCCPA <- c(MCCP[1:100], MCCP[401:500])
  MCCPB <- c(MCCP[101:200], MCCP[501:600])
  MCCPC <- c(MCCP[201:300], MCCP[601:700])
  MCCPD <- c(MCCP[301:400], MCCP[701:800])
  tCP <- c(trueCP[[1]], trueCP[[2]])
  if(nCP == 0) {
    range <- 1:100
  } else if(nCP == 1) {
    range <- 101:200
  } else {
    range <- 1:200
  }
  tOF <- as.table(rbind(score(MCCPA[range], tCP[range], radius, npts = 480),
                        score(MCCPB[range], tCP[range], radius, npts = 480),
                        score(MCCPC[range], tCP[range], radius, npts = 480),
                        score(MCCPD[range], tCP[range], radius, npts = 480),
                        score(lapply(NPMOJO[range], function(x) x[[1]]), tCP[range], radius, npts = 480),
                        score(lapply(NPMOJO[range], function(x) x[[2]]), tCP[range], radius, npts = 480),
                        score(lapply(NPMOJO[range], function(x) x[[3]]), tCP[range], radius, npts = 480),
                        score(lapply(CNPMOJO[range], function(x) x[[1]]), tCP[range], radius, npts = 480),
                        score(lapply(CNPMOJO[range], function(x) x[[2]]), tCP[range], radius, npts = 480),
                        score(lapply(CNPMOJO[range], function(x) x[[3]]), tCP[range], radius, npts = 480),
                        c(1, 0, 0, 0, 1)))
  dimnames(tOF) <- list(Method = c("MCCP[FWER]_q = 0.05", "MCCP[FWER]_q = 0.10",
                                   "MCCP[FDR]_q = 0.05", "MCCP[FDR]_q = 0.10",
                                   "NP-MOJO[G = 160]_q = 0.01", "NP-MOJO[G = 160]_q = 0.05", "NP-MOJO[G = 160]_q = 0.10",
                                   "CNP-MOJO[G = 40]_q = 0.01", "CNP-MOJO[G = 40]_q = 0.05", "CNP-MOJO[G = 40]_q = 0.10",
                                   paste0("Oracle_", process, cp)),
                        Metric = c("TPR", "FPR", "FDR", "FOR", "MCC"))
  return(tOF)
}

online_heatplot <- function(tableOBJ) {
  df <- data.frame(tableOBJ)
  df <- df %>% filter(Metric %in% c("TPR", "FDR", "MCC"))
  df$Freq[df$Metric == "FDR"] <- 1 - df$Freq[df$Metric == "FDR"]
  df$Metric <- factor(df$Metric, levels = c("MCC", "FDR", "PPV", "TPR"))
  df$Metric[df$Metric == "FDR"] <- "PPV"
  df$Parameters <- unlist(lapply(str_split(df$Method, "_"), function(x) x[2]))
  df$Method <- unlist(lapply(str_split(df$Method, "_"), function(x) x[1]))
  df$Metric <- factor(df$Metric, levels = c("MCC", "PPV", "TPR"))
  df$Method <- factor(df$Method, levels = c("Oracle", "SCCP[FWER]", "SCCP[FDR]", 
                                            "ScanB", "kCUSUM", 
                                            "MCCP[FWER]", "MCCP[FDR]", 
                                            "EDiv[MS = 30]", "EDiv[MS = 160]",
                                            "NP-MOJO[G = 160]", "KCP", "SBS1", "SBS2", 
                                            "FWER[q = 0.05]", "FWER[q = 0.10]", "FDR[q = 0.05]", "FDR[q = 0.10]"))
  df$Parameters <- factor(df$Parameters, levels = c("TAR", "TAR0", "TAR1", "TAR2", "GP", "GP0", "GP1", "GP2", "Naive", 
                                                    "q = 0.01", "q = 0.05", "q = 0.10",
                                                    "ARL = 1e+06", "ARL = 1e+05", "ARL = 1e+04", "ARL = 1e+03", "ARL = 1e+02", "ARL = 1e+01",
                                                    "C = 100", "C = 80", "C = 60", "C = 40", "C = 20", "C = 2",
                                                    "\u03b3 = 80", "\u03b3 = 161", "\u03b3 = 171", "\u03b3 = 181", "\u03b3 = 191"))
  df0 <- df[df$Method == "Oracle",]
  df1 <- df[df$Method %in% c("SCCP[FWER]", "ScanB"),]
  df2 <- df[df$Method %in% c("SCCP[FDR]", "kCUSUM"),]
  p0 <- ggplot() +
    geom_tile(data = df0, aes(x = Parameters, y = Metric, fill = Freq), color = "black") +
    facet_wrap(~ Method, nrow = 1, scales = "free_x", labeller = label_parsed) +
    scale_x_discrete("", position = "bottom", expand = c(0, 0)) +
    scale_y_discrete("", expand = c(0, 0)) +
    geom_text(data = df0, aes(x = Parameters, y = Metric, label = sprintf("%.2f", Freq), color = ifelse(Freq > 0.7, "white", "black")), size = 4) +
    scale_color_manual(values = c("white" = "white", "black" = "black"), guide = "none") +
    scale_fill_gradient2(low = "white", mid = "white", high = "gray30", midpoint = 1/5, limits = c(0, 1),
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", frame.linewidth = 0.5)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(angle = 45, hjust = 0.7),
          axis.text.x = element_text(angle = 45, hjust = 1, color = "white"),
          axis.title = element_blank(),
          axis.text = element_text(size = 12),
          legend.position = "right",
          legend.text = element_text(size = 14),
          legend.title = element_blank(),
          legend.key.height = unit(0.7, "in"),
          strip.background = element_blank(),
          strip.text = element_text(size = 16),
          panel.spacing = unit(1, "lines"),
          panel.border = element_rect(linewidth = 2))
  legend <- get_legend(p0)
  p0 <- p0 + theme(legend.position = "none")
  gp0 <- ggplotGrob(p0)
  facet.columns <- gp0$layout$l[grepl("panel", gp0$layout$name)]
  x.var0 <- sapply(ggplot_build(p0)$layout$panel_scales_x, function(l) length(l$range$range))
  gp0$widths[facet.columns] <- gp0$widths[facet.columns] * x.var0
  p1 <- ggplot() +
    geom_tile(data = df1, aes(x = Parameters, y = Metric, fill = Freq), color = "black") +
    facet_wrap(~ Method, nrow = 1, scales = "free_x", labeller = label_parsed) +
    scale_x_discrete("", position = "bottom", expand = c(0, 0)) +
    scale_y_discrete("", expand = c(0, 0)) +
    geom_text(data = df1, aes(x = Parameters, y = Metric, label = sprintf("%.2f", Freq), color = ifelse(Freq > 0.7, "white", "black")), size = 4) +
    scale_color_manual(values = c("white" = "white", "black" = "black"), guide = "none") +
    scale_fill_gradient2(low = "white", mid = "white", high = "gray30", midpoint = 1/5, limits = c(0, 1),
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", frame.linewidth = 0.5)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(angle = 45, hjust = 0.7),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank(),
          axis.text = element_text(size = 12),
          legend.position = "none",
          legend.key.width = unit(2, "cm"),
          strip.background = element_blank(),
          strip.text = element_text(size = 16),
          panel.spacing = unit(1, "lines"),
          panel.border = element_rect(linewidth = 2))
  gp1 <- ggplotGrob(p1)
  facet.columns <- gp1$layout$l[grepl("panel", gp1$layout$name)]
  x.var1 <- sapply(ggplot_build(p1)$layout$panel_scales_x, function(l) length(l$range$range))
  gp1$widths[facet.columns] <- gp1$widths[facet.columns] * x.var1
  p2 <- ggplot() +
    geom_tile(data = df2, aes(x = Parameters, y = Metric, fill = Freq), color = "black") +
    facet_wrap(~ Method, nrow = 1, scales = "free_x", labeller = label_parsed) +
    scale_x_discrete("", position = "bottom", expand = c(0, 0)) +
    scale_y_discrete("", expand = c(0, 0)) +
    geom_text(data = df2, aes(x = Parameters, y = Metric, label = sprintf("%.2f", Freq), color = ifelse(Freq > 0.7, "white", "black")), size = 4) +
    scale_color_manual(values = c("white" = "white", "black" = "black"), guide = "none") +
    scale_fill_gradient2(low = "white", mid = "white", high = "gray30", midpoint = 1/5, limits = c(0, 1),
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", frame.linewidth = 0.5)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(angle = 45, hjust = 0.7),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank(),
          axis.text = element_text(size = 12),
          legend.position = "none",
          legend.key.width = unit(2, "cm"),
          strip.background = element_blank(),
          strip.text = element_text(size = 16),
          panel.spacing = unit(1, "lines"),
          panel.border = element_rect(linewidth = 2))
  gp2 <- ggplotGrob(p2)
  facet.columns <- gp2$layout$l[grepl("panel", gp2$layout$name)]
  x.var2 <- sapply(ggplot_build(p2)$layout$panel_scales_x, function(l) length(l$range$range))
  gp2$widths[facet.columns] <- gp2$widths[facet.columns] * x.var2
  layout <- rbind(c(rep(NA, 5), rep(2, 20), rep(NA, 4)),
                  c(rep(1, 4), rep(NA, 1), rep(2, 20), rep(NA, 1), rep(4, 2), rep(NA, 1)),
                  c(rep(1, 4), rep(NA, 1), rep(3, 20), rep(NA, 1), rep(4, 2), rep(NA, 1)),
                  c(rep(NA, 5), rep(3, 20), rep(NA, 4)))
  gp <- arrangeGrob(gp0, gp1, gp2, legend, layout_matrix = layout)
  return(gp)
}

offline_heatplot <- function(tableOBJ) {
  df <- data.frame(tableOBJ)
  df <- df %>% filter(Metric %in% c("TPR", "FDR", "MCC"))
  df$Freq[df$Metric == "FDR"] <- 1 - df$Freq[df$Metric == "FDR"]
  df$Metric <- factor(df$Metric, levels = c("MCC", "FDR", "PPV", "TPR"))
  df$Metric[df$Metric == "FDR"] <- "PPV"
  df$Parameters <- unlist(lapply(str_split(df$Method, "_"), function(x) x[2]))
  df$Method <- unlist(lapply(str_split(df$Method, "_"), function(x) x[1]))
  df$Metric <- factor(df$Metric, levels = c("MCC", "PPV", "TPR"))
  df$Method <- factor(df$Method, levels = c("Oracle", "SCCP[FWER]", "SCCP[FDR]", 
                                            "ScanB", "kCUSUM", 
                                            "MCCP[FWER]", "MCCP[FDR]", 
                                            "EDiv[MS = 30]", "EDiv[MS = 160]",
                                            "NP-MOJO[G = 160]", "KCP", "SBS1", "SBS2", 
                                            "FWER[q = 0.05]", "FWER[q = 0.10]", "FDR[q = 0.05]", "FDR[q = 0.10]"))
  df$Parameters <- factor(df$Parameters, levels = c("TAR", "TAR0", "TAR1", "TAR2", "GP", "GP0", "GP1", "GP2", "Naive", 
                                                    "q = 0.01", "q = 0.05", "q = 0.10",
                                                    "ARL = 1e+06", "ARL = 1e+05", "ARL = 1e+04", "ARL = 1e+03", "ARL = 1e+02", "ARL = 1e+01",
                                                    "C = 100", "C = 80", "C = 60", "C = 40", "C = 20", "C = 2",
                                                    "\u03b3 = 80", "\u03b3 = 161", "\u03b3 = 171", "\u03b3 = 181", "\u03b3 = 191"))
  df0 <- df[df$Method == "Oracle",]
  df1 <- df[df$Method %in% c("MCCP[FWER]", "EDiv[MS = 160]", "SBS1", "SBS2"),]
  df2 <- df[df$Method %in% c("MCCP[FDR]", "NP-MOJO[G = 160]", "KCP"),]
  p0 <- ggplot() +
    geom_tile(data = df0, aes(x = Parameters, y = Metric, fill = Freq), color = "black") +
    facet_wrap(~ Method, nrow = 1, scales = "free_x", labeller = label_parsed) +
    scale_x_discrete("", position = "bottom", expand = c(0, 0)) +
    scale_y_discrete("", expand = c(0, 0)) +
    geom_text(data = df0, aes(x = Parameters, y = Metric, label = sprintf("%.2f", Freq), color = ifelse(Freq > 0.7, "white", "black")), size = 4) +
    scale_color_manual(values = c("white" = "white", "black" = "black"), guide = "none") +
    scale_fill_gradient2(low = "white", mid = "white", high = "gray30", midpoint = 1/5, limits = c(0, 1),
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", frame.linewidth = 0.5)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(angle = 45, hjust = 0.7),
          axis.text.x = element_text(angle = 45, hjust = 1, color = "white"),
          axis.title = element_blank(),
          axis.text = element_text(size = 12),
          legend.position = "right",
          legend.text = element_text(size = 14),
          legend.title = element_blank(),
          legend.key.height = unit(0.7, "in"),
          strip.background = element_blank(),
          strip.text = element_text(size = 16),
          panel.spacing = unit(1, "lines"),
          panel.border = element_rect(linewidth = 2))
  legend <- get_legend(p0)
  p0 <- p0 + theme(legend.position = "none")
  gp0 <- ggplotGrob(p0)
  facet.columns <- gp0$layout$l[grepl("panel", gp0$layout$name)]
  x.var0 <- sapply(ggplot_build(p0)$layout$panel_scales_x, function(l) length(l$range$range))
  gp0$widths[facet.columns] <- gp0$widths[facet.columns] * x.var0
  p1 <- ggplot() +
    geom_tile(data = df1, aes(x = Parameters, y = Metric, fill = Freq), color = "black") +
    facet_wrap(~ Method, nrow = 1, scales = "free_x", labeller = label_parsed) +
    scale_x_discrete("", position = "bottom", expand = c(0, 0)) +
    scale_y_discrete("", expand = c(0, 0)) +
    geom_text(data = df1, aes(x = Parameters, y = Metric, label = sprintf("%.2f", Freq), color = ifelse(Freq > 0.7, "white", "black")), size = 4) +
    scale_color_manual(values = c("white" = "white", "black" = "black"), guide = "none") +
    scale_fill_gradient2(low = "white", mid = "white", high = "gray30", midpoint = 1/5, limits = c(0, 1),
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", frame.linewidth = 0.5)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(angle = 45, hjust = 0.7),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank(),
          axis.text = element_text(size = 12),
          legend.position = "none",
          legend.key.width = unit(2, "cm"),
          strip.background = element_blank(),
          strip.text = element_text(size = 16),
          panel.spacing = unit(1, "lines"),
          panel.border = element_rect(linewidth = 2))
  gp1 <- ggplotGrob(p1)
  facet.columns <- gp1$layout$l[grepl("panel", gp1$layout$name)]
  x.var1 <- sapply(ggplot_build(p1)$layout$panel_scales_x, function(l) length(l$range$range))
  gp1$widths[facet.columns] <- gp1$widths[facet.columns] * x.var1
  p2 <- ggplot() +
    geom_tile(data = df2, aes(x = Parameters, y = Metric, fill = Freq), color = "black") +
    facet_wrap(~ Method, nrow = 1, scales = "free_x", labeller = label_parsed) +
    scale_x_discrete("", position = "bottom", expand = c(0, 0)) +
    scale_y_discrete("", expand = c(0, 0)) +
    geom_text(data = df2, aes(x = Parameters, y = Metric, label = sprintf("%.2f", Freq), color = ifelse(Freq > 0.7, "white", "black")), size = 4) +
    scale_color_manual(values = c("white" = "white", "black" = "black"), guide = "none") +
    scale_fill_gradient2(low = "white", mid = "white", high = "gray30", midpoint = 1/5, limits = c(0, 1),
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", frame.linewidth = 0.5)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(angle = 45, hjust = 0.7),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank(),
          axis.text = element_text(size = 12),
          legend.position = "none",
          legend.key.width = unit(2, "cm"),
          strip.background = element_blank(),
          strip.text = element_text(size = 16),
          panel.spacing = unit(1, "lines"),
          panel.border = element_rect(linewidth = 2))
  gp2 <- ggplotGrob(p2)
  facet.columns <- gp2$layout$l[grepl("panel", gp2$layout$name)]
  x.var2 <- sapply(ggplot_build(p2)$layout$panel_scales_x, function(l) length(l$range$range))
  gp2$widths[facet.columns] <- gp2$widths[facet.columns] * x.var2
  layout <- rbind(c(rep(NA, 5), rep(2, 28), rep(NA, 5)),
                  c(rep(1, 4), rep(NA, 1), rep(2, 28), rep(NA, 1), rep(4, 2), rep(NA, 2)),
                  c(rep(1, 4), rep(NA, 1), rep(3, 28), rep(NA, 1), rep(4, 2), rep(NA, 2)),
                  c(rep(NA, 5), rep(3, 28), rep(NA, 5)))
  gp <- arrangeGrob(gp0, gp1, gp2, legend, layout_matrix = layout)
  return(gp)
}

naive_heatplot <- function(tableOBJ) {
  df <- data.frame(tableOBJ)
  df <- df %>% filter(Metric %in% c("TPR", "FDR", "MCC"))
  df$Freq[df$Metric == "FDR"] <- 1 - df$Freq[df$Metric == "FDR"]
  df$Metric <- factor(df$Metric, levels = c("MCC", "FDR", "PPV", "TPR"))
  df$Metric[df$Metric == "FDR"] <- "PPV"
  df$Parameters <- unlist(lapply(str_split(df$Method, "_"), function(x) x[2]))
  df$Method <- unlist(lapply(str_split(df$Method, "_"), function(x) x[1]))
  df$Metric <- factor(df$Metric, levels = c("MCC", "PPV", "TPR"))
  df$Method <- factor(df$Method, levels = c("Oracle", "SCCP[FWER]", "SCCP[FDR]", 
                                            "ScanB", "kCUSUM", 
                                            "MCCP[FWER]", "MCCP[FDR]", 
                                            "EDiv[MS = 30]", "EDiv[MS = 160]",
                                            "NP-MOJO[G = 160]", "KCP", "SBS1", "SBS2", 
                                            "FWER[q = 0.05]", "FWER[q = 0.10]", "FDR[q = 0.05]", "FDR[q = 0.10]"))
  df$Parameters <- factor(df$Parameters, levels = c("TAR", "TAR0", "TAR1", "TAR2", "GP", "GP0", "GP1", "GP2", "Naive", 
                                                    "q = 0.01", "q = 0.05", "q = 0.10",
                                                    "ARL = 1e+06", "ARL = 1e+05", "ARL = 1e+04", "ARL = 1e+03", "ARL = 1e+02", "ARL = 1e+01",
                                                    "C = 100", "C = 80", "C = 60", "C = 40", "C = 20", "C = 2",
                                                    "\u03b3 = 80", "\u03b3 = 161", "\u03b3 = 171", "\u03b3 = 181", "\u03b3 = 191"))
  df0 <- df[df$Method == "Oracle",]
  df1 <- df[df$Method %in% c("FWER[q = 0.05]", "FWER[q = 0.10]"),]
  df2 <- df[df$Method %in% c("FDR[q = 0.05]", "FDR[q = 0.10]"),]
  p0 <- ggplot() +
    geom_tile(data = df0, aes(x = Parameters, y = Metric, fill = Freq), color = "black") +
    facet_wrap(~ Method, nrow = 1, scales = "free_x", labeller = label_parsed) +
    scale_x_discrete("", position = "bottom", expand = c(0, 0)) +
    scale_y_discrete("", expand = c(0, 0)) +
    geom_text(data = df0, aes(x = Parameters, y = Metric, label = sprintf("%.2f", Freq), color = ifelse(Freq > 0.7, "white", "black")), size = 4) +
    scale_color_manual(values = c("white" = "white", "black" = "black"), guide = "none") +
    scale_fill_gradient2(low = "white", mid = "white", high = "gray30", midpoint = 1/5, limits = c(0, 1),
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", frame.linewidth = 0.5)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(angle = 45, hjust = 0.7),
          axis.text.x = element_text(angle = 45, hjust = 1, color = "white"),
          axis.title = element_blank(),
          axis.text = element_text(size = 12),
          legend.position = "right",
          legend.text = element_text(size = 14),
          legend.title = element_blank(),
          legend.key.height = unit(0.7, "in"),
          strip.background = element_blank(),
          strip.text = element_text(size = 16),
          panel.spacing = unit(1, "lines"),
          panel.border = element_rect(linewidth = 2))
  legend <- get_legend(p0)
  p0 <- p0 + theme(legend.position = "none")
  gp0 <- ggplotGrob(p0)
  facet.columns <- gp0$layout$l[grepl("panel", gp0$layout$name)]
  x.var0 <- sapply(ggplot_build(p0)$layout$panel_scales_x, function(l) length(l$range$range))
  gp0$widths[facet.columns] <- gp0$widths[facet.columns] * x.var0
  p1 <- ggplot() +
    geom_tile(data = df1, aes(x = Parameters, y = Metric, fill = Freq), color = "black") +
    facet_wrap(~ Method, nrow = 1, scales = "free_x", labeller = label_parsed) +
    scale_x_discrete("", position = "bottom", expand = c(0, 0)) +
    scale_y_discrete("", expand = c(0, 0)) +
    geom_text(data = df1, aes(x = Parameters, y = Metric, label = sprintf("%.2f", Freq), color = ifelse(Freq > 0.7, "white", "black")), size = 4) +
    scale_color_manual(values = c("white" = "white", "black" = "black"), guide = "none") +
    scale_fill_gradient2(low = "white", mid = "white", high = "gray30", midpoint = 1/5, limits = c(0, 1),
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", frame.linewidth = 0.5)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(angle = 45, hjust = 0.7),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank(),
          axis.text = element_text(size = 12),
          legend.position = "none",
          legend.key.width = unit(2, "cm"),
          strip.background = element_blank(),
          strip.text = element_text(size = 16),
          panel.spacing = unit(1, "lines"),
          panel.border = element_rect(linewidth = 2))
  gp1 <- ggplotGrob(p1)
  facet.columns <- gp1$layout$l[grepl("panel", gp1$layout$name)]
  x.var1 <- sapply(ggplot_build(p1)$layout$panel_scales_x, function(l) length(l$range$range))
  gp1$widths[facet.columns] <- gp1$widths[facet.columns] * x.var1
  p2 <- ggplot() +
    geom_tile(data = df2, aes(x = Parameters, y = Metric, fill = Freq), color = "black") +
    facet_wrap(~ Method, nrow = 1, scales = "free_x", labeller = label_parsed) +
    scale_x_discrete("", position = "bottom", expand = c(0, 0)) +
    scale_y_discrete("", expand = c(0, 0)) +
    geom_text(data = df2, aes(x = Parameters, y = Metric, label = sprintf("%.2f", Freq), color = ifelse(Freq > 0.7, "white", "black")), size = 4) +
    scale_color_manual(values = c("white" = "white", "black" = "black"), guide = "none") +
    scale_fill_gradient2(low = "white", mid = "white", high = "gray30", midpoint = 1/5, limits = c(0, 1),
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", frame.linewidth = 0.5)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(angle = 45, hjust = 0.7),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank(),
          axis.text = element_text(size = 12),
          legend.position = "none",
          legend.key.width = unit(2, "cm"),
          strip.background = element_blank(),
          strip.text = element_text(size = 16),
          panel.spacing = unit(1, "lines"),
          panel.border = element_rect(linewidth = 2))
  gp2 <- ggplotGrob(p2)
  facet.columns <- gp2$layout$l[grepl("panel", gp2$layout$name)]
  x.var2 <- sapply(ggplot_build(p2)$layout$panel_scales_x, function(l) length(l$range$range))
  gp2$widths[facet.columns] <- gp2$widths[facet.columns] * x.var2
  layout <- rbind(c(rep(NA, 5), rep(2, 27), rep(NA, 5)),
                  c(rep(1, 4), rep(NA, 1), rep(2, 27), rep(NA, 1), rep(4, 2), rep(NA, 2)),
                  c(rep(1, 4), rep(NA, 1), rep(3, 27), rep(NA, 1), rep(4, 2), rep(NA, 2)),
                  c(rep(NA, 5), rep(3, 27), rep(NA, 5)))
  gp <- arrangeGrob(gp0, gp1, gp2, legend, layout_matrix = layout)
  return(gp)
}

Cheatplot <- function(tableOBJ) {
  df <- data.frame(tableOBJ)
  df <- df %>% filter(Metric %in% c("TPR", "FDR", "MCC"))
  df$Freq[df$Metric == "FDR"] <- 1 - df$Freq[df$Metric == "FDR"]
  df$Metric <- factor(df$Metric, levels = c("MCC", "FDR", "PPV", "TPR"))
  df$Metric[df$Metric == "FDR"] <- "PPV"
  df$Parameters <- unlist(lapply(str_split(df$Method, "_"), function(x) x[2]))
  df$Method <- unlist(lapply(str_split(df$Method, "_"), function(x) x[1]))
  df$Metric <- factor(df$Metric, levels = c("MCC", "PPV", "TPR"))
  df$Method <- factor(df$Method, levels = c("Oracle",
                                            "MCCP[FWER]", "MCCP[FDR]",
                                            "NP-MOJO[G = 160]", "CNP-MOJO[G = 40]"))
  df$Parameters <- factor(df$Parameters, levels = c("", "q = 0.01", "q = 0.05", "q = 0.10"))
  df$Parameters[is.na(df$Parameters)] <- ""
  p0 <- ggplot() +
    geom_tile(data = df, aes(x = Parameters, y = Metric, fill = Freq), color = "black") +
    facet_wrap(~ Method, nrow = 1, scales = "free_x", labeller = label_parsed) +
    scale_x_discrete("", position = "bottom", expand = c(0, 0)) +
    scale_y_discrete("", expand = c(0, 0)) +
    geom_text(data = df, aes(x = Parameters, y = Metric, label = sprintf("%.2f", Freq), color = ifelse(Freq > 0.7, "white", "black")), size = 4) +
    scale_color_manual(values = c("white" = "white", "black" = "black"), guide = "none") +
    scale_fill_gradient2(low = "white", mid = "white", high = "gray30", midpoint = 1/5, limits = c(0, 1),
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", frame.linewidth = 0.5)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(angle = 45, hjust = 0.7),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank(),
          axis.text = element_text(size = 12),
          legend.position = "right",
          legend.text = element_text(size = 14),
          legend.title = element_blank(),
          legend.key.height = unit(0.7, "in"),
          strip.background = element_blank(),
          strip.text = element_text(size = 16),
          panel.spacing = unit(1, "lines"),
          panel.border = element_rect(linewidth = 2))
  legend <- get_legend(p0)
  p0 <- p0 + theme(legend.position = "none")
  gp0 <- ggplotGrob(p0)
  facet.columns <- gp0$layout$l[grepl("panel", gp0$layout$name)]
  x.var0 <- sapply(ggplot_build(p0)$layout$panel_scales_x, function(l) length(l$range$range))
  gp0$widths[facet.columns] <- gp0$widths[facet.columns] * x.var0
  return(gp0)
}


