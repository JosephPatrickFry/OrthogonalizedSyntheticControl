#Estimate tuning parameter for control weights
EstimateLambdaDelta <- function(Y0, YJ, Z, T1) {
  
  Y0 = as.matrix(Y0)
  YJ = as.matrix(YJ)
  Z = as.matrix(Z)
  
  T0 = ncol(as.matrix(YJ))
  J = nrow(as.matrix(YJ))
  K = nrow(as.matrix(Z))
  
  H <- diag(J) 
  C <- rep(0,J)
  A <- rbind(Z %*% t(YJ)/sqrt(T0), -Z %*% t(YJ)/sqrt(T0))
  b <- c(Z %*% t(Y0)/sqrt(T0), -Z %*% t(Y0)/sqrt(T0))
  
  E = t(as.matrix(c(0,rep(1,J))))
  f = as.vector(1)
  G = cbind(as.matrix(rep(1,2*K)),A)
  h = b
  
  res = limSolve::linp(E = E, F = f, G = G, H = h, Cost = c(1,rep(0,J)),
       verbose = FALSE,lower = 0, upper = NULL)
  
  lambdaDelta = res$X[1]
  lambdaDelta = lambdaDelta*log(min(c(T0,T1)))*log(J)/log(K)
  
  return(lambdaDelta)
}

#Estimate tuning parameter for instrument weights
EstimateLambdaEta <- function(PreY0, PreYJ, PostY0, PostYJ, Z) {
  
  PreY0 = t(as.matrix(PreY0))
  PreYJ = as.matrix(PreYJ)
  Z = as.matrix(Z)
  
  PostY0 = t(as.matrix(PostY0))
  PostYJ = as.matrix(PostYJ)
  
  T0 = ncol(as.matrix(PreYJ))
  T1 = nrow(as.matrix(PostYJ))

  gdelta = rbind(Z %*% t(PreYJ)/T0, colMeans(PostYJ))
  
  K = nrow(gdelta)
  J = ncol(gdelta)

  E = t(as.matrix(c(rep(0,K),1)))
  f = as.vector(1)
  G = cbind(as.matrix(rep(1,2*J)),rbind(t(gdelta), - t(gdelta)))
  h = c(rep(0, J*2))
  
  res = limSolve::linp(E = E, F = f, G = G, H = h, Cost = c(1,rep(0,K)),
             verbose = FALSE,lower = NULL, upper = NULL, ispos = FALSE)
  lambdaEta = res$X[1]*log(min(c(T0,T1)))
  return(lambdaEta)
}

#Estimate control weights
EstimateDelta <- function(Y0, YJ, Z, Scaled = FALSE) {
  #PreY0: treated unit's pretreatment outcomes
  #PreYJ: control units' pretreatment outcomes
  #Z: instrument units' pretreatment outcomes
  #Scaled: Rescales outcomes by time periods if true

  Y0 = t(as.matrix(Y0))
  YJ = t(as.matrix(YJ))
  Z = t(as.matrix(Z))
  
  T0 = nrow(as.matrix(YJ))
  J = ncol(as.matrix(YJ))
  K = ncol(as.matrix(Z))

  if (Scaled) {
    big.dataframe = cbind(Y0, YJ, Z)
    divisor <- sqrt(apply(big.dataframe, 1, var))
    data.scaled <- t(t(big.dataframe) %*% (1/(divisor) * diag(rep(dim(big.dataframe)[1],1))))
  
    
    Y0.scaled <- t(as.matrix(data.scaled[, c(1:(dim(Y0)[2]))]))
    YJ.scaled <- t(data.scaled[, c((dim(Y0)[2]+1):(dim(Y0)[2]+dim(YJ)[2]))])
    Z.scaled <- t(as.matrix(data.scaled[, c((dim(Y0)[2]+dim(YJ)[2]+1):(dim(Y0)[2]+dim(YJ)[2]+dim(Z)[2]))]))
    
  } else {
    Y0.scaled = t(as.matrix(Y0))
    YJ.scaled = t(as.matrix(YJ))
    Z.scaled = t(as.matrix(Z))
  }
  

  lambda = EstimateLambdaDelta(Y0 = Y0.scaled, YJ = YJ.scaled, Z = Z.scaled, T1 = T1)
  
  A <- rbind(Z.scaled %*% t(YJ.scaled)/sqrt(T0), -Z.scaled %*% t(YJ.scaled)/sqrt(T0))
  h <- c(Z.scaled %*% t(Y0.scaled)/sqrt(T0) - lambda, -Z.scaled %*% t(Y0.scaled)/sqrt(T0) - lambda)
  
  res = limSolve::ldei(E = t(as.matrix(rep(1,J))), F = as.vector(1), G = A, H = h,lower = 0, upper = 1,verbose = FALSE)
  delta = res$X

  return(list(delta = delta, lambda = lambda))
  
}

#Estimate instrument weights
EstimateNormalizedEta <- function(PreY0, PreYJ, PostY0, PostYJ, Z, Scaled = FALSE) {
  #PreY0: treated unit's pretreatment outcomes
  #PreYJ: control units' pretreatment outcomes
  #Z: instrument units' pretreatment outcomes
  #PostY0: treated unit's posttreatment outcomes
  #PostYJ: control units' posttreatment outcomes
  #Scaled: Rescales outcomes by time periods if true
  
  PreY0 = t(as.matrix(PreY0))
  PreYJ = t(as.matrix(PreYJ))
  Z = t(as.matrix(Z))
  
  PostY0 = t(as.matrix(PostY0))
  PostYJ = t(as.matrix(PostYJ))
  
  T0 = nrow(as.matrix(PreY0))
  T1 = nrow(as.matrix(PostY0))
  
  if (Scaled) {
    big.dataframe = cbind(PreY0, PreYJ, Z)
    divisor <- sqrt(apply(big.dataframe, 1, var))
    data.scaled <- t(t(big.dataframe) %*% (1/(divisor) * diag(rep(dim(big.dataframe)[1],1))))
    PreY0.scaled <- t(as.matrix(data.scaled[, c(1:(dim(PreY0)[2]))]))
    PreYJ.scaled <- t(data.scaled[, c((dim(PreY0)[2]+1):(dim(PreY0)[2]+dim(PreYJ)[2]))])
    Z.scaled <- t(as.matrix(data.scaled[, c((dim(PreY0)[2]+dim(PreYJ)[2]+1):(dim(PreY0)[2]+dim(PreYJ)[2]+dim(Z)[2]))]))
    
    big.dataframe = cbind(PostY0, PostYJ)
    divisor <- sqrt(apply(big.dataframe, 1, var))
    data.scaled <- t(t(big.dataframe) %*% (1/(divisor) * diag(rep(dim(big.dataframe)[1],1))))
    PostY0.scaled <- as.matrix(data.scaled[, c(1:(dim(PostY0)[2]))])
    PostYJ.scaled <- data.scaled[, c((dim(PostY0)[2]+1):(dim(PostY0)[2]+dim(PostYJ)[2]))]
    
  } else {
    PreY0.scaled = t(as.matrix(PreY0))
    PreYJ.scaled = t(as.matrix(PreYJ))
    Z.scaled = t(as.matrix(Z))
    
    PostY0.scaled = as.matrix(PostY0)
    PostYJ.scaled = as.matrix(PostYJ)
  }
  
  gdelta = rbind(Z.scaled %*% t(PreYJ.scaled)/T0, colMeans(PostYJ.scaled))
  
  K = nrow(gdelta)
  J = ncol(gdelta)
  initial = c(rep(0,K-1),1)
  A = rbind(t(gdelta), - t(gdelta))

  lambda = EstimateLambdaEta(PreY0 = PreY0.scaled, PreYJ = PreYJ.scaled, 
                               Z = Z.scaled, PostY0 = PostY0.scaled, PostYJ = PostYJ.scaled)
  h = c(rep(-lambda, J*2))
  res = limSolve::ldei(E = t(as.matrix(c(rep(0,K-1),1))), F = as.vector(1),G = A, H = h,
             lower = NA, upper = NA,verbose = FALSE)

  eta = res$X
  
  return(list(eta = eta, lambda = lambda))
  
}
