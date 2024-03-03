library(tidyverse)
library(progress)
library(matrixStats)
library(gridExtra)
library(cowplot)
library(ecp)
library(MASS)
library(hdbinseg)
library(abind)
library(Rcpp)
library(RcppArmadillo)

genchanges <- function(length, n_changes, start_point, end_point, min_spacing) {
  if(n_changes > (end_point - start_point) / min_spacing) {
    n_changes = floor((end_point - start_point) / min_spacing)
  }
  if(n_changes > 1) {
    change_flag <- FALSE
    while(!change_flag) {
      change_points <- sort(sample(start_point:end_point, n_changes, replace = F))
      if(min(diff(change_points)) >= min_spacing) {
        change_flag <- TRUE
      }
    }
  } else if(n_changes == 1) {
    change_points <- sort(sample(start_point:end_point, n_changes, replace = F))
  } else {
    change_points <- NULL
  }
  change_points <- c(0, change_points, length)
  return(change_points)
}

genTARmat <- function(unif) {
  tol <- 1
  A <- array(dim = c(4, 4, 4))
  reg <- 1
  count <- 0
  while(reg <= 4) {
    err <- 1
    a111 <- ifelse(reg != 1, runif(1, -unif, unif), 0)
    a211 <- ifelse(reg != 1, runif(1, -unif, unif), 0)
    a121 <- ifelse(reg != 3, runif(1, -unif, unif), 0)
    a221 <- ifelse(reg != 3, runif(1, -unif, unif), 0)
    a112 <- ifelse(reg != 4, runif(1, -unif, unif), 0)
    a212 <- ifelse(reg != 4, runif(1, -unif, unif), 0)
    a122 <- ifelse(reg != 2, runif(1, -unif, unif), 0)
    a222 <- ifelse(reg != 2, runif(1, -unif, unif), 0)
    A[,,reg] <- round(rbind(matrix(c(a111, a121, a112, a122, a211, a221, a212, a222), nrow = 2), cbind(diag(2), matrix(0, nrow = 2, ncol = 2))), 2)
    A[,,reg][abs(A[,,reg]) < 0.2] <- sign(A[,,reg][abs(A[,,reg]) < 0.2]) * 0.2
    evs <- abs(eigen(A[,,1])$values)
    err <- max(abs(abs(eigen(A[,,reg])$values) - evs))
    count <- count + 1
    if(err < 0.01) {
      reg <- reg + 1
    }
    if((max(evs) > 0.9) | (min(evs) < 0.05)) {
      reg <- 1
      count <- 0
    }
    if(count > 10000) {
      reg <- 1
      count <- 0
    }
  }
  return(A)
}

genGP <- function(length = 800, n_changes = 0, min_spacing) {
  start_point <- min_spacing
  end_point <- length - min_spacing + 1
  change_points <- genchanges(length, n_changes, start_point, end_point, min_spacing)
  freq <- runif(1, 10, 30)
  SM <- matrix(nrow = length, ncol = length)
  for(t1 in seq(length)) {
    for(t2 in seq(length)) {
      SM[t1, t2] = exp(-32 * (sin((t1 - t2) / freq) ^ 2))
      SM[t2, t1] = SM[t1, t2]
    }
  }
  corr <- runif(1, -0.8, 0.8)
  SIG <- kronecker(SM, matrix(c(1, corr, corr, 1), nrow = 2))
  out0 <- mvrnorm(mu = rep(0, 2 * length), Sigma = SIG)
  out1 <- cbind(out0[seq(1, 2 * length, 2)], out0[seq(2, 2 * length, 2)])
  out2 <- vector("list", length(change_points) - 1)
  for(section in 1:(length(change_points) - 1)) {
    out2[[section]] <- (-1) ^ (section - 1) * out1[(change_points[section] + 1):change_points[section + 1], ]
  }
  signal <- do.call(rbind, lapply(out2, scale))
  return(list(data = signal, change_points = change_points[(change_points > 0) & (change_points < length)], freq = freq, corr = corr))
}

genTAR <- function(length = 800, n_changes = 0, min_spacing, TARmat) {
  start_point <- min_spacing
  end_point <- length - min_spacing + 1
  change_points <- genchanges(length, n_changes, start_point, end_point, min_spacing)
  state <- vector("numeric", length = length)
  regime <- vector("numeric", length = length)
  output <- matrix(0, nrow = length + 1, 8)
  error <- cbind(matrix(rnorm(2 * length), ncol = 2), matrix(0, nrow = length, ncol = 2))
  output[1,] <- rnorm(8)
  for(time in 1:length) {
    state[time] <- (sum(time > change_points) - 1) %% 2
    regime[time] <- ifelse(sum(output[time,]) < 0, 0, 1)
    output[time + 1, 1:4] <- t(TARmat[,, 2 * (state[time] + 1) - (1 - regime[time])] %*% matrix(output[time, 1:4]) + error[time,])
    output[time + 1, 5:8] <- output[time, 3:6]
  }
  output <- output[2:(length + 1), 1:2]
  return(list(data = output, change_points = change_points[(change_points > 0) & (change_points < length)], state = state, regime = regime, TAR = TARmat))
}

genNaive <- function(length = 800, spacing) {
  if(spacing == 80) {
    change_points <- c(240, 320, 400, 480, 560)
  } else if(spacing == 161) {
    change_points <- c(220, 381, 542)
  } else if(spacing == 171) {
    change_points <- c(220, 391, 562)
  } else if(spacing == 181) {
    change_points <- c(220, 401, 582)
  } else if(spacing == 191) {
    change_points <- c(220, 411, 602)
  }
  cps <- c(change_points, length)
  out <- matrix(rnorm(2 * change_points[1]), ncol = 2)
  for(i in 2:length(cps)) {
    out <- rbind(out, matrix(rnorm(2 * (cps[i] - cps[i - 1]), sd = (1 - 0.8 * ((i + 1) %% 2))), ncol = 2))
  }
  return(list(data = out, change_points = change_points))
}
