Simulations = function(i = NULL, J = 14, K = 9, T1 = 16, T0_values = c(30,60,90), loadings = loadings,
                       Factor_pars = Factor_pars, Sigma = Sigma, lambda_delta = 0, beta0 = 0) {
  source("FM_DGP.R")
  source("SeriesHAC.R")
  source("GMMSC.R")
  source("SubsamplingInferenceSC.R")
  library(MSCMT)
  library(rngtools)
  library(mvtnorm)
  library(Synth)
  library(data.table)
  library(synthdid)
  library(comprehenr)
  library(doParallel)
  library(mvtnorm)
  library(gmm)
  library(geometry)
  library(forecast)
  library(matlib)
  library(scinference)
  set.seed(i)
  N = ncol(loadings)
  T0_total = T0_values[length(T0_values)]
  data = FM_DGP(N = N,Time = T0_total+T1,loadings = loadings,
                Factor_pars = Factor_pars, Sigma = Sigma) 
  data$assigment = (data$year > T0_total) &(data$Countryno <= 5)
  sim = panel.matrices(data, unit = 1, time = 2, outcome = 3, treatment = 5)
  Y = t(sim$Y) 
  
  is_control = c(TRUE,TRUE,TRUE,TRUE,FALSE,TRUE,FALSE,TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,FALSE,TRUE,TRUE,TRUE,FALSE,TRUE,FALSE,TRUE)
  is_treated = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE)
  is_instrument = c(FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,TRUE,FALSE,FALSE,TRUE,TRUE,FALSE,TRUE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE)

  TotalY0 = as.matrix(Y[,is_treated])
  TotalYJ = as.matrix(Y[,is_control])
  TotalZ = as.matrix(Y[,is_instrument]) 
  
  PostYJ = t(as.matrix(TotalYJ[(T0_total+1):(T0_total+T1),]))
  PostY0 = t(as.matrix(TotalY0[(T0_total+1):(T0_total+T1),])) + beta0
  
  Results = rep(NA, 90)
  
  
  for (t in 0:(length(T0_values)-1)) {
    T0 = T0_values[t+1]
    PreYJ = t(as.matrix(TotalYJ[(T0_total+1-T0):T0_total,]))
    Z = t(as.matrix(TotalZ[(T0_total+1-T0):T0_total,]))
    Z = rbind(Z,rep(1,T0))
    PreY0 = t(as.matrix(TotalY0[(T0_total+1-T0):T0_total,]))
    
    
    #My Method
    res = EstimateDelta(Y0 = PreY0, YJ = PreYJ, Z = Z, constrained = TRUE, lambda = lambda_delta, quadopt = "Idei")
    delta = res$delta
    #lambda_eta = log(J)*res$lambda
    
    #eta = EstimateEta(PreY0 = PreY0, PreYJ = PreYJ, PostY0 = PostY0, PostYJ = PostYJ, Z = Z,constrained = TRUE, 
    #                  delta = delta, lambda = lambda_eta)
    res = EstimateNormalizedEta(PreY0 = PreY0, PreYJ = PreYJ, PostY0 = PostY0, PostYJ = PostYJ, Z = Z,constrained = TRUE, 
                                delta = delta, lambda = lambda_delta)
    eta = res$eta 
    lambda_eta = res$lambda
    g0 = as.matrix(rbind(Z %*% t(PreY0)/T0, mean(PostY0)))
    #print(eta)
    #print(g0)
    gdelta = rbind(Z %*% t(PreYJ)/T0, rowMeans(PostYJ))
    beta = as.numeric(eta %*% (g0 - gdelta %*% delta)/eta[K])
    #print(beta)
    Postg = PostY0 - t(t(PostYJ) %*% delta) - beta
    Preg = matrix(NA, nrow = nrow(Z), ncol = T0)
    for (s in 1:T0) {
      Preg[,s] = t(Z[,s])*as.numeric((PreY0[s] - t(PreYJ[,s]) %*% delta))
    }
    #hOptimal = 2
    hOptimal = Optimalh(Preg = Preg, Postg = Postg, p = 1, alpha = alpha)
    #V = SeriesHAC(Preg, Postg, eta, h = hOptimal, alpha = alpha)
    V = SeriesHACFast(Preg, Postg, eta, h = hOptimal)
    #V = 1.5
    pvalue = pvaluetTest(beta0 = 0, betahat = beta, V = V,h = hOptimal,n = min(c(T0,T1)))
    
    Results[1+30*t] = (beta -beta0)^2 #MSE
    Results[2+30*t] = beta - beta0 #Bias
    Results[3+30*t] = pvalue #pvalue
    Results[4+30*t] = hOptimal #Amount of Smoothing
    Results[5+30*t] = V #Variance
    
    
    wOLS = GMMSC(Y0 = PreY0, YJ = PreYJ, YK = diag(T0), Efficient = FALSE, quadopt = "LowRankQP", meanfit = FALSE)$Weights
    betaOLS = mean(t(PostY0) - t(PostYJ) %*% wOLS)
    A_m = ResamplingSC_OLS(m = floor(T0/3), iter = 100, W = wOLS, Y0 = t(PreY0), YJ = t(PreYJ), PostY0 = t(PostY0), PostYJ = t(PostYJ),
                             Replace = FALSE, LongRunVar = FALSE)
  
    pvalueOLS = pValueSC(H0 = 0, stat = betaOLS, stat_subsamples = A_m)
    Results[6+30*t] = betaOLS^2 #MSE
    Results[7+30*t] = betaOLS #Bias
    Results[8+30*t] = pvalueOLS #pvalue
    #Results[9+30*t] = norm(wOLS,type = "2") #WNorm
    #Results[10+30*t] = mean((PostY0 - PostYJ %*% wOLS))^2 #ATT_MSE
    
    Kfold = 3
    res = scinference(Y1 = as.vector(c(PreY0,PostY0)), Y0 = t(cbind(PreYJ,PostYJ)), T1 = T1, T0 = T0,
                      inference_method = "ttest",alpha = 0.1, theta0 = 0,
                      estimation_method = "sc",permutation_method = "mb",K = Kfold) 
    tstat = sqrt(Kfold)*(res$att)/res$se  
    ttestpvalue = 2*pt(-abs(tstat), df = Kfold-1)
    Results[11+30*t] = (res$att - beta0)^2 #MSE
    Results[12+30*t] = res$att - beta0 #Bias
    Results[13+30*t] = ttestpvalue #pvalue
    #Results[14+30*t] = norm(wOLS,type = "2") #WNorm
    #Results[15+30*t] = mean((PostY0 - PostYJ %*% wOLS))^2 #ATT_MSE
    
    #SDID
    #Ysdid = cbind(rbind(PreYJ,PreY0),rbind(PostYJ,PostY0))
    #tau.hat = synthdid_estimate(Y = Ysdid, N0 = J, T0 = T0)
    #se = sqrt(vcov(tau.hat, method='placebo'))
    #SDIDpvalue = 2*pnorm(-abs(tau.hat/se))
    
    #Results[26+30*t] = tau.hat^2 #MSE
    #Results[27+30*t] = tau.hat #Bias
    #Results[28+30*t] = SDIDpvalue #pvalue
    #Results[29+30*t] = norm(wOLS,type = "2") #WNorm
    #Results[30+30*t] = mean((PostY0 - PostYJ %*% wOLS))^2 #ATT_MSE
  }
  
  return(Results)
}

SimulationsInference = function(i = NULL, J = 14, K = 9, T1 = 16, T0_values = c(30,60,90), loadings = loadings,
                       Factor_pars = Factor_pars, Sigma = Sigma, lambda_delta = 0, beta0 = 0) {
  source("FM_DGP.R")
  source("SeriesHAC.R")
  source("GMMSC.R")
  source("SubsamplingInferenceSC.R")
  library(MSCMT)
  library(rngtools)
  library(mvtnorm)
  library(Synth)
  library(data.table)
  library(synthdid)
  library(comprehenr)
  library(doParallel)
  library(mvtnorm)
  library(gmm)
  library(geometry)
  library(forecast)
  library(matlib)
  library(scinference)
  set.seed(i)
  N = ncol(loadings)
  T0_total = T0_values[length(T0_values)]
  data = FM_DGP(N = N,Time = T0_total+T1,loadings = loadings,
                Factor_pars = Factor_pars, Sigma = Sigma) 
  data$assigment = (data$year > T0_total) &(data$Countryno <= 5)
  sim = panel.matrices(data, unit = 1, time = 2, outcome = 3, treatment = 5)
  Y = t(sim$Y) 
  
  is_control = c(TRUE,TRUE,TRUE,TRUE,FALSE,TRUE,FALSE,TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,FALSE,TRUE,TRUE,TRUE,FALSE,TRUE,FALSE,TRUE)
  is_treated = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE)
  is_instrument = c(FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,TRUE,FALSE,FALSE,TRUE,TRUE,FALSE,TRUE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE)
  
  TotalY0 = as.matrix(Y[,is_treated])
  TotalYJ = as.matrix(Y[,is_control])
  TotalZ = as.matrix(Y[,is_instrument]) 
  
  PostYJ = t(as.matrix(TotalYJ[(T0_total+1):(T0_total+T1),]))
  PostY0 = t(as.matrix(TotalY0[(T0_total+1):(T0_total+T1),]))
  
  Results = rep(NA, 90)
  
  
  for (t in 0:(length(T0_values)-1)) {
    T0 = T0_values[t+1]
    PreYJ = t(as.matrix(TotalYJ[(T0_total+1-T0):T0_total,]))
    Z = t(as.matrix(TotalZ[(T0_total+1-T0):T0_total,]))
    Z = rbind(Z,rep(1,T0))
    PreY0 = t(as.matrix(TotalY0[(T0_total+1-T0):T0_total,]))
    
    
    #My Method
    res = EstimateDelta(Y0 = PreY0, YJ = PreYJ, Z = Z, constrained = TRUE, lambda = lambda_delta, quadopt = "Idei")
    delta = res$delta
    #lambda_eta = log(J)*res$lambda
    
    #eta = EstimateEta(PreY0 = PreY0, PreYJ = PreYJ, PostY0 = PostY0, PostYJ = PostYJ, Z = Z,constrained = TRUE, 
    #                  delta = delta, lambda = lambda_eta)
    res = EstimateNormalizedEta(PreY0 = PreY0, PreYJ = PreYJ, PostY0 = PostY0, PostYJ = PostYJ, Z = Z,constrained = TRUE, 
                                delta = delta, lambda = lambda_delta)
    eta = res$eta 
    lambda_eta = res$lambda
    g0 = as.matrix(rbind(Z %*% t(PreY0)/T0, mean(PostY0)))
    gdelta = rbind(Z %*% t(PreYJ)/T0, rowMeans(PostYJ))
    beta = as.numeric(eta %*% (g0 - gdelta %*% delta)/eta[K])
    Postg = PostY0 - t(t(PostYJ) %*% delta) - beta
    Preg = matrix(NA, nrow = nrow(Z), ncol = T0)
    for (s in 1:T0) {
      Preg[,s] = t(Z[,s])*as.numeric((PreY0[s] - t(PreYJ[,s]) %*% delta))
    }
    hOptimal = Optimalh(Preg = Preg, Postg = Postg, p = 1, alpha = alpha)
    V = SeriesHACFast(Preg, Postg, eta, h = hOptimal)
    #V = SeriesHAC(Preg, Postg, eta, h = hOptimal, alpha = alpha)
    pvalue = pvaluetTest(beta0 = beta0, betahat = beta, V = V,h = hOptimal,n = min(c(T0,T1)))
    
    Results[1+30*t] = beta^2 #MSE
    Results[2+30*t] = beta #Bias
    Results[3+5*t] = pvalue #pvalue
    Results[4+5*t] = hOptimal #Amount of Smoothing
    Results[5+30*t] = V #ATT-MSE
    
    
    #wOLS = GMMSC(Y0 = PreY0, YJ = PreYJ, YK = diag(T0), Efficient = FALSE, quadopt = "LowRankQP", meanfit = FALSE)$Weights
    #betaOLS = mean(t(PostY0) - t(PostYJ) %*% wOLS)
    #A_m = ResamplingSC_OLS(m = floor(T0/3), iter = 100, W = wOLS, Y0 = t(PreY0), YJ = t(PreYJ), PostY0 = t(PostY0), PostYJ = t(PostYJ),
    #                         Replace = TRUE, LongRunVar = TRUE)
    
    #pvalueOLS = pValueSC(H0 = beta0, stat = betaOLS, stat_subsamples = A_m)
    #Results[6+30*t] = betaOLS^2 #MSE
    #Results[7+30*t] = betaOLS #Bias
    #Results[8+30*t] = pvalueOLS #pvalue
    #Results[9+30*t] = norm(wOLS,type = "2") #WNorm
    #Results[10+30*t] = mean((PostY0 - PostYJ %*% wOLS))^2 #ATT_MSE
    
    #t-test
    Kfold = 3
    res = scinference(Y1 = as.vector(c(PreY0,PostY0)), Y0 = t(cbind(PreYJ,PostYJ)), T1 = T1, T0 = T0,
                      inference_method = "ttest",alpha = 0.1, theta0 = beta0,
                      estimation_method = "sc",permutation_method = "mb",K = Kfold) 
    tstat = sqrt(Kfold)*(res$att)/res$se  
    ttestpvalue = 2*pt(-abs(tstat), df = Kfold-1)
    Results[11+30*t] = res$att^2 #MSE
    Results[12+30*t] = res$att #Bias
    Results[13+30*t] = ttestpvalue #pvalue
    #Results[14+30*t] = norm(wOLS,type = "2") #WNorm
    #Results[15+30*t] = mean((PostY0 - PostYJ %*% wOLS))^2 #ATT_MSE
    
    #Conformal Inference
    res = scinference(Y1 = as.vector(c(PreY0,PostY0)), Y0 = t(cbind(PreYJ,PostYJ)), T1 = T1, T0 = T0,
                      inference_method = "conformal",alpha = 0.1,theta0 = beta0,
                      estimation_method = "sc",permutation_method = "mb")
    Conpvalue = res$p_val
    
    Results[16+30*t] = res$att^2 #MSE
    Results[17+30*t] = res$att #Bias
    Results[18+30*t] = Conpvalue #pvalue
    #Results[19+30*t] = norm(wOLS,type = "2") #WNorm
    #Results[20+30*t] = mean((PostY0 - PostYJ %*% wOLS))^2 #ATT_MSE
    
    #End-of-the-Sample Instability Test
    source("CaoDowd2019functions.R")
    res = sp_andrews_te(Y0 = t(rbind(PreY0,PreYJ)),Y1 = rbind(PostY0,PostYJ),
                        A = as.matrix(c(1,rep(0,J))),loo=F,Normal=TRUE,alpha_sig=0.05)
    
    Results[21+30*t] = mean(res$Estimates)^2 #MSE
    Results[22+30*t] = mean(res$Estimates) #Bias
    Results[23+30*t] = res$pval #pvalue
    #Results[24+30*t] = norm(wOLS,type = "2") #WNorm
    #Results[25+30*t] = mean((PostY0 - PostYJ %*% wOLS))^2 #ATT_MSE
    
    #SDID
    Ysdid = cbind(rbind(PreYJ,PreY0),rbind(PostYJ,PostY0))
    tau.hat = synthdid_estimate(Y = Ysdid, N0 = J, T0 = T0)
    se = sqrt(vcov(tau.hat, method='placebo'))
    SDIDpvalue = 2*pnorm(-abs(tau.hat/se))
    
    Results[26+30*t] = tau.hat^2 #MSE
    Results[27+30*t] = tau.hat #Bias
    Results[28+30*t] = SDIDpvalue #pvalue
    #Results[29+30*t] = norm(wOLS,type = "2") #WNorm
    #Results[30+30*t] = mean((PostY0 - PostYJ %*% wOLS))^2 #ATT_MSE
  }
  
  return(Results)
}


