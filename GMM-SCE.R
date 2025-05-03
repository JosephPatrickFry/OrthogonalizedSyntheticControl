#Function to Estimate the Control Weights with One-Step GMM
GMMSC <- function(Y0 = NULL, YJ = NULL, YK = NULL, meanfit = TRUE) {
  #Y0: treated unit's pretreatment outcomes
  #YJ: control units' pretreatment outcomes
  #YK: instrument units' pretreatment outcomes
  #meanfit: Addtionally uses a constant as an instrument if true
  
  
  #Loading packages
  library(kernlab)
  library(LowRankQP)
  library(matlib)
  library(cointReg)
  library(quadprog)
  library(corpcor)
  library(sandwich)
  
  Y0 = as.matrix(Y0)
  YJ = as.matrix(YJ)
  YK = as.matrix(YK)
  
  T0 = nrow(as.matrix(Y0))
  J = ncol(as.matrix(YJ))
  if(meanfit == TRUE) {
    YK <- cbind(rep(1,T0), YK)
  }
  K = ncol(as.matrix(YK))
  
  #Normalize data
  big.dataframe = cbind(Y0, YJ, YK)
  divisor <- sqrt(apply(big.dataframe, 1, var))
  data.scaled <- t(t(big.dataframe) %*% (1/(divisor) * diag(rep(dim(big.dataframe)[1],1))))
  Y0.scaled <- as.matrix(data.scaled[, c(1:(dim(Y0)[2]))])
  YJ.scaled <- data.scaled[, c((dim(Y0)[2]+1):(dim(Y0)[2]+dim(YJ)[2]))]
  YK.scaled <- as.matrix(data.scaled[, c((dim(Y0)[2]+dim(YJ)[2]+1):(dim(Y0)[2]+dim(YJ)[2]+dim(YK)[2]))])
  
  
  H <- t(YJ.scaled) %*% YK.scaled %*% t(YK.scaled) %*% YJ.scaled 
  c <- as.matrix(-1*c(t(Y0.scaled) %*% YK.scaled %*% t(YK.scaled) %*% YJ.scaled))
  A <- t(rep(1, J))
  b <- 1
  l <- rep(0, J)
  u <- rep(1, J)
  r <- 0
  res <- LowRankQP(Vmat = H, dvec = c, Amat = A, bvec = 1, uvec = rep(1, length(c)), method = "LU")
  solution.w <- as.matrix(res$alpha)
  
  

  #Calculate J-statistic
  g = t(YK.scaled) %*% (Y0.scaled - YJ.scaled %*% solution.w)
  Jstat <- as.numeric(t(g) %*% diag(K) %*% g *1/T0)
  
  
  #Return control weights and J-statistic
  return(list(Weights = solution.w, Jstatistic = Jstat))
}


#Function to Select the Set of Control Units
ModelSelection <- function(Y0 = NULL, YN0 = NULL, YN1 = NULL, pvalue = .01) {
  #Y0: treated unit pretreatment outcomes 
  #YN0: pretreatment outcomes of units which may either be controls or instruments
  #YN1: pretreatment outcomes of units which are guaranteed to be used as instruments
  
  library(comprehenr)
  YN0 = as.matrix(YN0)
  N0 = ncol(as.matrix(YN0))
  T0 = nrow(as.matrix(YN0))
  
  #Initialize the number of controls
  extraJ = 1
  
  #Calculate pvalue for Sargan–Hansen test
  pvalue = exp(-sqrt(T0))
  
  
  #Preform Initial Estimate using all units in YN0 as instruments
  if (is.null(YN1)) {
    result = GMMSC(Y0 = Y0, YJ = YN0, YK = YN0)
    N1 = 0
  } else {
    result = GMMSC(Y0 = Y0, YJ = YN0, YK = cbind(YN0,YN1))
    N1 = ncol(YN1)
  }
  testStat = result$Jstatistic
  InitialControls = result$Weights
  
  #Obtain initial set of controls based on units which recieve non-zero weight
  JInitial = (InitialControls > .00001)
  YJinitial = YN0[,JInitial]
  YN0 = YN0[,!JInitial]
  J = 1 + ncol(YJinitial)
  criticalValue = qchisq(pvalue, df=max(N1+1+N0 -2*J,1), lower.tail = FALSE)
  rejected = (testStat > criticalValue)
  if (!rejected) {
    return(JInitial)
  }


  rejected = TRUE
  
  #Sort units based on mean-squared difference to treated unit
  if (N0 > 1) {
    mse = to_vec(for(i in 1:N0) mean((Y0 - YN0[,i])^2))
    sorted = sort(mse, decreasing = FALSE)
  }
  
  while ((rejected == TRUE) & (extraJ < N0 + as.numeric(!is.null(YN1)))) {
    if (N0 > 1) {
      controls = mse %in% sorted[1:extraJ]
    } else {
      controls = c(TRUE)
    }
    if (is.null(YJinitial)) {
      YJ = as.matrix(YN0[,controls])
    } else {
      YJ = cbind(YJinitial, YN0[,controls])
    }
    if (is.null(YN1)) {
      YK = as.matrix(YN0[,!controls])
    } else {
      YK = cbind(YN1, YN0[,!controls])
    }
    
    K = ncol(YK)+1
    
    testStat = GMMSC(Y0 = Y0, YJ = YJ, YK = YK)$Jstatistic
    
    #Critical value from chi-squared distribution
    criticalValue = qchisq(pvalue, df=max(K-J,1), lower.tail = FALSE)

    #If rejecting the test, switch one unit from instrument to control and repeat
    #Else return this set of of controls
    rejected = (testStat > criticalValue)
    extraJ = extraJ + 1
    J = J + 1
    K = K - 1
  }
  return(controls)
}
