source("~/Desktop/Research/Rodu/MCCP/utils3.R")
sourceCpp("~/Desktop/Research/Rodu/MCCP/MCCP.cpp")
load("~/Desktop/Research/Rodu/MCCP/Data/TARMat.RData")

cnp.mojo.ESM <- function(data, washL, trainL, threshold) {
  L <- nrow(data)
  data_dim <- ncol(data)
  params <- RNNParamFit(data, washL, trainL, 5)
  W1 <- genW1(params[1], data_dim, 100)
  RS0 <- RNN1(data, W1, 0, params[2], params[3])
  C <- CCalc(RS0[,, (washL + 1):(washL + trainL)], 100)
  CRS <- CRNN1(data, W1, C, 0, washL, params[2], params[3])
  angles <- angleCalc(CRS[,, 1:L], CRS[,, (L + 1):(2 * L)])
  seq <- apply(angles, 1, mean)[(washL + trainL + 1):L]
  ests <- np.mojo.multilag(seq, G = washL, lags = 0:2, alpha = threshold)$cpts[,1] + washL + trainL
  if(length(ests) == 0) {
    est <- L
  } else {
    est <- ests
  }
  return(est)
}

cnp.mojo.SCCP <- function(data, washL, trainL, threshold, FWER) {
  L <- nrow(data)
  data_dim <- ncol(data)
  Pmax <- floor((L - washL - trainL) / (washL + trainL + 1))
  if(FWER) {
    cutoff <- threshold / seq(Pmax, 1, -1)
  } else {
    cutoff <- threshold * rep(1, Pmax)
  }
  n_rejects <- 0
  fwd <- cnp.mojo.ESM(data, washL, trainL, cutoff[1])
  ests <- fwd
  while(fwd <= (L - 2 * washL - 2 * trainL + 1)) {
    n_rejects <- n_rejects + 1
    newdata <- data[(fwd + 1):L,]
    fwd <- cnp.mojo.ESM(newdata, washL, trainL, cutoff[n_rejects + 1]) + fwd
    ests <- c(ests, fwd)
  }
  return(ests)
}

cnp.mojo.MCCP <- function(data, washL, trainL, threshold, FWER) {
  L <- nrow(data)
  fwd <- cnp.mojo.SCCP(data, washL, trainL, threshold / 2, FWER) 
  bwd <- L - cnp.mojo.SCCP(data[nrow(data):1,], washL, trainL, threshold / 2, FWER)
  tau <- Reconcile(matrix(fwd), matrix(bwd), washL + trainL + 1, L)
  tau <- tau[(tau > 0) & (tau < L)]
  return(tau)
}