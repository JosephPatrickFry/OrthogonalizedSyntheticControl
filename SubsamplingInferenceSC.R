BlockResamplingSC_OLS <- function(m, iter, FUN = NULL, Y0, YJ, PostY0, PostYJ, quadopt = "LowRankQP", 
                              Efficient = FALSE, LongRunVar = TRUE) {
  library(kernlab)
  library(LowRankQP)
  library(matlib)
  library(cointReg)
  library(stats)
  
  T0 = nrow(YJ)
  T1 = nrow(PostYJ)
  J = ncol(PostYJ)
  
  synth.out = synth(X0 = YJ, X1 = Y0, Z0 = YJ, Z1 = Y0, custom.v = rep(1/(T0), T0), quadopt = quadopt)
  W = synth.out$solution.w

  v_hat = PostY0 - PostYJ %*% W
  if (LongRunVar) {
    Sigma_v = getLongRunVar(v_hat)$Omega
  } else {
    Sigma_v = var(v_hat)
  }
  W_M = c()
  for (i in 1:(T0 - m + 1)) {
    n = i:(i+m-1)
    Y0_m = as.matrix(Y0[n,])
    YJ_m = YJ[n,]
    synth.out_m = synth(X0 = YJ_m, X1 = Y0_m, Z0 = YJ_m, Z1 = Y0_m, custom.v = rep(1/(m), m), quadopt = quadopt)
    W_m = synth.out_m$solution.w 
    W_M = c(W_M, W_m)
  }
  result = c()
  for (i in 1: iter) {
    s = sample(T0-m+1, 1) - 1
    W_m = W_M[(s*J + 1):(s*J + J)]
    #v_m = rnorm(n = T1, mean = 0, sd = sqrt(Sigma_v))
    v_m = rt(n = 1, df = T1-1)*sqrt(Sigma_v)
    A_m = -sqrt(T1/T0)*t(colMeans(PostYJ)) %*% (sqrt(m)*(W_m - W)) + v_m#sum(v_m)/sqrt(T1)
    result = c(result, A_m)
  }
  
  
  return(result)
}

ResamplingSC_OLS <- function(m, iter, Y0, YJ, PostY0, PostYJ, quadopt = "LowRankQP", 
                         W = NULL, LongRunVar = FALSE, Replace = TRUE) {
  library(kernlab)
  library(LowRankQP)
  library(matlib)
  library(cointReg)
  library(stats)
  source("GMMSC.R")
  
  T0 = nrow(YJ)
  T1 = nrow(PostYJ)
  
  #synth.out = synth(X0 = YJ, X1 = Y0, Z0 = YJ, Z1 = Y0, custom.v = rep(1/(T0), T0), quadopt = quadopt)
  #W = synth.out$solution.w
  
  v_hat = PostY0 - PostYJ %*% W
  if (LongRunVar) {
    Sigma_v = getLongRunVar(v_hat)$Omega
  } else {
    Sigma_v = var(v_hat)
  }
  
  result = c()
  for (i in 1:iter) {
    n = sample(T0, m, replace = Replace)
    Y0_m = as.matrix(Y0[n,])
    YJ_m = YJ[n,]
    W_m = GMMSC(Y0 = t(Y0_m), YJ = t(YJ_m), YK = diag(m), Efficient = FALSE, quadopt = "LowRankQP", meanfit = FALSE)$Weights
    #synth.out_m = synth(X0 = YJ_m, X1 = Y0_m, Z0 = YJ_m, Z1 = Y0_m, custom.v = rep(1/(m), m), quadopt = quadopt)
    #W_m = synth.out_m$solution.w    
    #v_m = rnorm(n = T1, mean = 0, sd = sqrt(Sigma_v))
    v_m = rt(n = 1, df = T1-1)*sqrt(Sigma_v)
    A_m = -sqrt(T1/T0)*t(colMeans(PostYJ)) %*% (sqrt(m)*(W_m - W)) + v_m#sum(v_m)/sqrt(T1)
    result = c(result, A_m)
  }
  return(result)
}


ResamplingSC <- function(m, iter, Y0, YJ, PostY0, PostYJ, YK, quadopt = "LowRankQP", 
                  Efficient = FALSE, LongRunVar = FALSE, Replace = TRUE) {
  library(kernlab)
  library(LowRankQP)
  library(matlib)
  library(cointReg)
  library(stats)
  
  T0 = nrow(YJ)
  T1 = nrow(PostYJ)
  
  W = GMMSC(Y0 = Y0, YJ = YJ, YK = YK, Efficient = Efficient, quadopt = quadopt)
  v_hat = PostY0 - PostYJ %*% W
  if (LongRunVar) {
    Sigma_v = getLongRunVar(v_hat)$Omega
  } else {
    Sigma_v = var(v_hat)
  }
  
  result = c()
  for (i in 1:iter) {
    n = sample(T0, m, replace = Replace)
    Y0_m = as.matrix(Y0[n,])
    YJ_m = YJ[n,]
    YK_m = YK[n,]
    W_m =  GMMSC(Y0 = Y0_m, YJ = YJ_m, YK = YK_m, Efficient = Efficient, quadopt = quadopt)
    #v_m = rnorm(n = T1, mean = 0, sd = sqrt(Sigma_v))
    v_m = rt(n = 1, df = T1-1)*sqrt(Sigma_v)
    A_m = -sqrt(T1/T0)*t(colMeans(PostYJ)) %*% (sqrt(m)*(W_m - W)) + v_m#sum(v_m)/sqrt(T1)
    result = c(result, A_m)
  }
  return(result)
}

ConfidenceIntervalsSC <- function(stat, stat_subsamples, alpha = .05, T1 = 100) {
  LB = stat - quantile(stat_subsamples, 1 - alpha/2)/sqrt(T1)
  UB = stat - quantile(stat_subsamples, alpha/2)/sqrt(T1)
  return(c(LB, UB))
}

pValueSC <- function(H0 = 0, stat, stat_subsamples) {
  dist = abs(H0 - stat_subsamples)
  diff = abs(H0 - stat)
  pvalue = mean(dist > diff)
  return(pvalue)
}

BlockResamplingSC <- function(m, iter, FUN = NULL, Y0, YJ, PostY0, PostYJ, YK, quadopt = "LowRankQP", 
                         Efficient = FALSE, LongRunVar = TRUE) {
  library(kernlab)
  library(LowRankQP)
  library(matlib)
  library(cointReg)
  library(stats)
  
  T0 = nrow(YJ)
  T1 = nrow(PostYJ)
  J = ncol(PostYJ)

  W = GMMSC(Y0 = Y0, YJ = YJ, YK = YK, Efficient = Efficient, quadopt = quadopt)
  v_hat = PostY0 - PostYJ %*% W
  if (LongRunVar) {
    Sigma_v = getLongRunVar(v_hat)$Omega
  } else {
    Sigma_v = var(v_hat)
  }
  W_M = c()
  for (i in 1:(T0 - m + 1)) {
    n = i:(i+m-1)
    Y0_m = as.matrix(Y0[n,])
    YJ_m = YJ[n,]
    YK_m = YK[n,]
    W_m =  GMMSC(Y0 = Y0_m, YJ = YJ_m, YK = YK_m, Efficient = Efficient, quadopt = quadopt)
    W_M = c(W_M, W_m)
  }
  result = c()
  for (i in 1: iter) {
    s = sample(T0-m+1, 1) - 1
    W_m = W_M[(s*J + 1):(s*J + J)]
    #v_m = rnorm(n = T1, mean = 0, sd = sqrt(Sigma_v))
    v_m = rt(n = 1, df = T1-1)*sqrt(Sigma_v)
    #t = sample(T1-m+1, 1) 
    #v_m = v_hat[t:(t+m-1)]
    #v_m = sample(v_hat, T1, replace = FALSE)
    A_m = -sqrt(T1/T0)*t(colMeans(PostYJ)) %*% (sqrt(m)*(W_m - W)) + v_m#sum(v_m)/sqrt(T1)
    result = c(result, A_m)
  }
  return(result)
}

DoubleBlockResamplingSC <- function(m, b, iter, FUN = NULL, Y0, YJ, PostY0, PostYJ, YK, quadopt = "LowRankQP", 
                              Efficient = FALSE, LongRunVar = TRUE) {
  library(kernlab)
  library(LowRankQP)
  library(matlib)
  library(cointReg)
  library(stats)
  
  T0 = nrow(YJ)
  T1 = nrow(PostYJ)
  J = ncol(PostYJ)
  W = GMMSC(Y0 = Y0, YJ = YJ, YK = YK, Efficient = Efficient, quadopt = quadopt)
  
  W_M = c()
  for (i in 1:(T0 - m + 1)) {
    n = i:(i+m-1)
    Y0_m = as.matrix(Y0[n,])
    YJ_m = YJ[n,]
    YK_m = YK[n,]
    W_m =  GMMSC(Y0 = Y0_m, YJ = YJ_m, YK = YK_m, Efficient = Efficient, quadopt = quadopt)
    W_M = c(W_M, W_m)
  }
  
  result = c()
  for (i in 1: iter) {
    s = sample(T0-m+1, 1) - 1
    W_m = W_M[(s*J + 1):(s*J + J)]
    t = sample(T1-b+1, 1)
    PostY0_b = as.matrix(PostY0[t:(t+b-1),])
    PostYJ_b = PostYJ[t:(t+b-1),]
    v_m = PostY0_b - PostYJ_b %*% W_m
    #v_m = sample(v_hat, T1, replace = FALSE)
    A_m = -sqrt(T1/T0)*t(colMeans(PostYJ)) %*% (sqrt(m)*(W_m - W)) + sum(v_m)/sqrt(b)
    result = c(result, A_m)
  }
  return(result)
}


ProjectionInference <- function(iter, Y0, YJ, PostY0, PostYJ, YK, quadopt = "LowRankQP", 
                         Efficient = FALSE, LongRunVar = TRUE) {
  library(kernlab)
  library(LowRankQP)
  library(matlib)
  library(cointReg)
  library(stats)
  library(Matrix)
  
  T0 = nrow(YJ)
  T1 = nrow(PostYJ)
  J = ncol(YJ)
  
  W = GMMSC(Y0 = Y0, YJ = YJ, YK = YK, Efficient = Efficient, quadopt = quadopt)
  v_hat = PostY0 - PostYJ %*% W
  if (LongRunVar) {
    Sigma_v = getLongRunVar(v_hat)$Omega
  } else {
    Sigma_v = mean((v_hat)^2)
  }
  W_0 = GMMSC(Y0 = Y0, YJ = YJ, YK = YK, Efficient = FALSE, quadopt = quadopt)
  G <- NULL
  Moments = cbind(rep(1,T0), YK)
  for (i in 1:ncol(Moments)) {
    G <- cbind(G, Moments[,i]*(Y0 - YJ %*% W_0))
  }
  if (Efficient) {
    Omega = getLongRunVar(G)$Omega
    divisor = norm(Omega, type = "F")/nrow(Omega)
    Omega = Omega/divisor
    test <<- nearPD(Omega, ensureSymmetry = TRUE, doSym = TRUE, conv.tol = 1e-10, base.matrix = TRUE, posd.tol = 1e-20)
    Omega = test$mat
    A_T <- inv(Omega)/divisor
  } else {
    A_T <- inv(getLongRunVar(G)$Sigma)
  }
  Variance <- t(YJ) %*% Moments %*% A_T %*% t(Moments) %*% YJ / T0^2
  divisor <- norm(Variance, type = "F")/nrow(Variance)
  Variance = inv(Variance)
  Variance = nearPD(Variance, ensureSymmetry = TRUE, base.matrix = TRUE)$mat
  H = t(YJ) %*% Moments %*% A_T %*% t(Moments) %*% YJ / (T0^2)
  result = c()
  for (i in 1:iter) {
    Z_T = rmvnorm(n=1, mean = rep(0,J), sigma = Variance)
    c = -1/(T0^2) * c(Z_T %*% t(YJ) %*% Moments %*% A_T %*% t(Moments) %*% YJ)
    Z_m =  GMMProjection(Z_T = Z_T, H = H, c = c)
    v_m = rnorm(n = T1, mean = 0, sd = sqrt(Sigma_v))
    A_m = -sqrt(T1/T0)*t(colMeans(PostYJ)) %*% (Z_m) + sum(v_m)/sqrt(T1)
    result = c(result, A_m)
  }
  return(result)
}
