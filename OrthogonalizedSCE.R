#Estimate of the average treatment effect on the treated unit which is 
#orthogonalized with respect to the control weights
OrthoganilzedSCE <- function(PreY0, PreYJ, Z, PostY0, PostYJ, alpha = .05, 
                             beta0 = 0, includeConstant = TRUE) {
  source("RegularizedEstimate.R")
  source("SeriesHAC.R")
  #PreY0: treated unit's pretreatment outcomes
  #PreYJ: control units' pretreatment outcomes
  #Z: instrument units' pretreatment outcomes
  #PostY0: treated unit's posttreatment outcomes
  #PostYJ: control units' posttreatment outcomes
  #includeConstant: Addtionally uses a constant as an instrument if true
  #beta0: Null value of average treatment effect on the treated unit
  #alpha: Confidence interval will be calculated for the 1-alpha level
  
  
  if(includeConstant) {
    Z = rbind(Z, rep(1,T0))
  }
  T0 = ncol(PreYJ)
  J = nrow(PreYJ)
  T1 = ncol(PostYJ)
  Q = nrow(Z)+1
  
  #Estimate Control Weights
  res = EstimateDelta(Y0 = PreY0, YJ = PreYJ, Z = Z, Scaled = TRUE)
  delta = res$delta
  
  #Estimate Instrument Weights
  res = EstimateNormalizedEta(PreY0 = PreY0, PreYJ = PreYJ, PostY0 = PostY0,
                              PostYJ = PostYJ, Z = Z, Scaled = TRUE)
  eta = res$eta 

  #Orthogonalized Moment Conditions
  g0 = as.matrix(rbind(Z %*% t(PreY0)/T0, mean(PostY0)))
  gdelta = rbind(Z %*% t(PreYJ)/T0, rowMeans(PostYJ))
  
  #Orthogonalized ATT
  beta = as.numeric(eta %*% (g0 - gdelta %*% delta)/eta[Q] )

  Postg = PostY0 - t(t(PostYJ) %*% delta) - beta
  Preg = matrix(NA, nrow = Q-1, ncol = T0)
  for (t in 1:T0) {
    Preg[,t] = t(Z[,t])*as.numeric((PreY0[t] - t(PreYJ[,t]) %*% delta))
  }
  #Calculate smoothing parameter using method of Sun (2013)
  hOptimal = CPEOptimalh(Preg = Preg, Postg = Postg, p = 1, sig = .05)
  #Estimate variance
  V = SeriesHACFast(Preg, Postg, eta, h = hOptimal)
  
  pvalue = as.numeric(pvaluetTest(beta0 = beta0, betahat = beta, V = V,
                                  h = hOptimal,n = min(c(T0,T1))))
  CI = CI_tTest(betahat = beta, V = V, h = hOptimal, alpha = alpha)
  
  return(list(beta = beta, pvalue = pvalue, CI = CI, df = hOptimal, ControlWeights = delta, InstrumentWeights = eta))
}
