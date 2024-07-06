EstimateDelta <- function(Y0 = NULL, YJ = NULL, Z = NULL, quadopt = "solve.QP", constrained = TRUE, 
                  Efficient = FALSE, lambda = 0) {
  library(kernlab)
  library(matlib)
  library(cointReg)
  library(quadprog)
  library(corpcor)
  library(sandwich)
  library(limSolve)
  
  Y0 = t(as.matrix(Y0))
  YJ = t(as.matrix(YJ))
  Z = t(as.matrix(Z))
  
  T0 = nrow(as.matrix(YJ))
  J = ncol(as.matrix(YJ))
  K = ncol(as.matrix(Z))

  big.dataframe = cbind(Y0, YJ, Z)
  divisor <- sqrt(apply(big.dataframe, 1, var))
  data.scaled <- t(t(big.dataframe) %*% (1/(divisor) * diag(rep(dim(big.dataframe)[1],1))))
  Y0.scaled <- t(as.matrix(data.scaled[, c(1:(dim(Y0)[2]))]))
  YJ.scaled <- t(data.scaled[, c((dim(Y0)[2]+1):(dim(Y0)[2]+dim(YJ)[2]))])
  Z.scaled <- t(as.matrix(data.scaled[, c((dim(Y0)[2]+dim(YJ)[2]+1):(dim(Y0)[2]+dim(YJ)[2]+dim(Z)[2]))]))
  
  H <- diag(J) 
  C <- rep(0,J)
  A <- rbind(Z.scaled %*% t(YJ.scaled)/sqrt(T0), -Z.scaled %*% t(YJ.scaled)/sqrt(T0))
  b <- c(Z.scaled %*% t(Y0.scaled)/sqrt(T0) - lambda, -Z.scaled %*% t(Y0.scaled)/sqrt(T0) - lambda)

  if (constrained) {   
    if (quadopt == "Idei") {
      #b <- c(Z.scaled %*% t(Y0.scaled)/sqrt(T0), -Z.scaled %*% t(Y0.scaled)/sqrt(T0))
      #lambda = max(abs(b - A %*% rep(1/J, J)))/10
      #print(lambda)
      #b <- c(Z.scaled %*% t(Y0.scaled)/sqrt(T0) - lambda, -Z.scaled %*% t(Y0.scaled)/sqrt(T0) - lambda)
      count = 0
      tempdelta = as.matrix(rep(1/J,J))
      while((round(norm(tempdelta,type = "2"),3) == round(sqrt(1/J),3))&count<30) {
        #print("test1")
        res = ldei(E = t(as.matrix(rep(1,J))), F = as.vector(1), G = A, H = b,lower = 0, upper = 1,verbose = FALSE)
        tempdelta = as.matrix(res$X)
        lambda = lambda/2
        b <- c(Z.scaled %*% t(Y0.scaled)/sqrt(T0) - lambda, -Z.scaled %*% t(Y0.scaled)/sqrt(T0) - lambda)
        count= count+1
      }
      count = 0
      while((round(norm(tempdelta,type = "2"),3) != round(sqrt(1/J),3))&count<30) {
        #print("test2")
        res = ldei(E = t(as.matrix(rep(1,J))), F = as.vector(1), G = A, H = b,lower = 0, upper = 1,verbose = FALSE)
        tempdelta = as.matrix(res$X)
        lambda = lambda/2
        b <- c(Z.scaled %*% t(Y0.scaled)/sqrt(T0) - lambda, -Z.scaled %*% t(Y0.scaled)/sqrt(T0) - lambda)
        count =count+1
      }
      lambda = lambda*4
      b <- c(Z.scaled %*% t(Y0.scaled)/sqrt(T0) - lambda, -Z.scaled %*% t(Y0.scaled)/sqrt(T0) - lambda)
      res = ldei(E = t(as.matrix(rep(1,J))), F = as.vector(1), G = A, H = b,lower = 0, upper = 1,verbose = FALSE)
      delta = res$X
      
    }
    b <- c(b,rep(0,J))
    A <- rbind(A, diag(J))
    b <- c(1,b)
    A <- rbind(t(as.matrix(rep(1,J))),A)
    if (quadopt == "solve.QP") {
      res <- solve.QP(Dmat = H, dvec = C, Amat = t(A), 
                      bvec = b, meq = 1)
      delta <- as.matrix(res$solution)
    }
  } else {
    if (quadopt == "solve.QP") {
      res <- solve.QP(Dmat = H, dvec = C, Amat = t(A), 
                      bvec = b, meq = 0)
      delta <- as.matrix(res$solution)
    }
  }

  return(list(delta = delta, lambda = lambda))
  
}

EstimateNormalizedEta <- function(PreY0, PreYJ, PostY0, PostYJ, Z, quadopt = "solve.QP", constrained = TRUE, 
                        delta = NULL, lambda = 0) {
  library(kernlab)
  library(matlib)
  library(cointReg)
  library(quadprog)
  library(corpcor)
  library(sandwich)
  
  PreY0 = t(as.matrix(PreY0))
  PreYJ = t(as.matrix(PreYJ))
  Z = t(as.matrix(Z))
  
  PostY0 = t(as.matrix(PostY0))
  PostYJ = t(as.matrix(PostYJ))
  
  T0 = nrow(as.matrix(PreY0))
  T1 = nrow(as.matrix(PostY0))
  J = ncol(as.matrix(PreYJ))
  K = ncol(as.matrix(Z))
  
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
  
  #beta = mean(PostY0.scaled - PostYJ.scaled %*% as.matrix(delta))
  
  #g0 = as.matrix(rbind(Z.scaled %*% PreY0.scaled/sqrt(T0), (mean(PostY0.scaled)-beta)*sqrt(T1)))
  gdelta = rbind(Z.scaled %*% t(PreYJ.scaled)/T0, colMeans(PostYJ.scaled))
  
  #initial = pseudoinverse(t(gdelta) %*% gdelta) %*% t(gdelta) %*% g0
  K = nrow(gdelta)
  J = ncol(gdelta)
  initial = c(rep(0,K-1),1)
  
  A = rbind(t(gdelta), - t(gdelta))
  #b = rep(colMeans(PostYJ.scaled)*sqrt(T1)-lambda,J)
  #upper = colMeans(PostYJ.scaled)
  #lower = -colMeans(PostYJ.scaled)
  #b = c(upper,lower) - lambda
  b = c(rep(-lambda, J*2))
  tempeta = rep(1,K)
  count = 0
  while((tempeta[1] != 0) & count<30) {
    #print("test1")
    res = ldei(E = t(as.matrix(c(rep(0,K-1),1))), F = as.vector(1),G = A, H = b,lower = NA, upper = NA,verbose = FALSE)
    tempeta = as.matrix(res$X)
    lambda = lambda/2
    b = c(rep(-lambda, J*2))
    count= count+1
  }
  lambda = lambda*4
  b = c(rep(-lambda, J*2))
  res = ldei(E = t(as.matrix(c(rep(0,K-1),1))), F = as.vector(1),G = A, H = b,lower = NA, upper = NA,verbose = FALSE)
  
  #res = solve.QP(Dmat = diag(K), dvec = rep(0,K), Amat = t(A), bvec = b)
  #res = ldei(E = t(as.matrix(rep(1,J))), F = as.vector(1), G = A, H = b,lower = 0, upper = 1,verbose = FALSE)
  
  #output = constrOptim(theta = initial, f = fr_eta, grad = grr_eta, ui = ui,
  #                     ci = ci, method = "BFGS")
  #eta <- optim(par = rep(1/K, K),fn = fr_eta, gr = grr_eta, method = "BFGS")$par
  
  
  #eta = output$par
  #eta = res$solution
  eta = res$X
  
  return(list(eta = eta, lambda = lambda))
  
}



EstimateEta <- function(PreY0, PreYJ, PostY0, PostYJ, Z, quadopt = "solve.QP", constrained = TRUE, 
                          delta = NULL, lambda = 0) {
  library(kernlab)
  library(matlib)
  library(cointReg)
  library(quadprog)
  library(corpcor)
  library(sandwich)
  
  PreY0 = t(as.matrix(PreY0))
  PreYJ = t(as.matrix(PreYJ))
  Z = t(as.matrix(Z))
  
  PostY0 = t(as.matrix(PostY0))
  PostYJ = t(as.matrix(PostYJ))
  
  T0 = nrow(as.matrix(PreY0))
  T1 = nrow(as.matrix(PostY0))
  J = ncol(as.matrix(PreYJ))
  K = ncol(as.matrix(Z))
  
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

  #beta = mean(PostY0.scaled - PostYJ.scaled %*% as.matrix(delta))
  
  #g0 = as.matrix(rbind(Z.scaled %*% PreY0.scaled/sqrt(T0), (mean(PostY0.scaled)-beta)*sqrt(T1)))
  gdelta = rbind(Z.scaled %*% t(PreYJ.scaled)/sqrt(T0), colMeans(PostYJ.scaled)*sqrt(T1))
  
  #initial = pseudoinverse(t(gdelta) %*% gdelta) %*% t(gdelta) %*% g0
  K = nrow(gdelta)
  eta_K = lambda/(100*max(abs(gdelta[,K])))
  initial = c(rep(0,K-1),eta_K)
  
  ui = rbind(t(gdelta), - t(gdelta))
  ci = c(rep(-lambda,J),rep(-lambda,J))
  

  output = constrOptim(theta = initial, f = fr_eta, grad = grr_eta, ui = ui,
                       ci = ci, method = "BFGS")
  #eta <- optim(par = rep(1/K, K),fn = fr_eta, gr = grr_eta, method = "BFGS")$par
  
  
  eta = output$par

  
  return(eta)
  
}

fr_eta <- function(x) {
  eta = x
  K = length(eta)
  Vg = matrix(1, nrow = K, ncol = K)
  as.numeric(t(as.matrix(eta)) %*% Vg %*% as.matrix(eta)/eta[K]^2)
  #as.matrix(eta) %*% Vg(beta,delta) %*% t(as.matrix(eta))/eta[K]^2
}

grr_eta <- function(x) { ## Gradient of 'fr'
  eta = x
  K = length(eta)
  Vg = matrix(1, nrow = K, ncol = K)
  
  temp = 2*sum(eta)/eta[K]^2
  
  Hi = as.numeric(t(eta) %*% Vg %*% eta)
  dHi = 2*sum(eta)
  Low = eta[K]^2
  dLow = 2*eta[K]
  temp2 = (Low*dHi - Hi*dLow)/Low^2
  
  c(rep(temp,K-1),temp2)
}

