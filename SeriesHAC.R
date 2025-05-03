phi = function(x,j) {
  if (j%%2 == 0) {
    return(sqrt(2)*sinpi(2*j*x))
  } else {
    return(sqrt(2)*cospi(2*j*x))
  }
}

Q = function(i,s,h) {
  temp = 0
  for(j in 0:(h-1)) {
    temp = temp + phi(i,j)*phi(s,j)
  }
  return(temp/h)
}

pvaluetTest = function(beta0 = 0, betahat, V, h,n) {
  t_n = sqrt(n)*(betahat - beta0)/sqrt(V)
  pvalue = 2*pt(-abs(t_n), df = h)
  return(pvalue)
}

CI_tTest <- function(betahat, V, h, alpha) {
  cv = abs(qt(alpha, df = h))
  CI = c(betahat - sqrt(V)*cv, betahat + sqrt(V)*cv)
  return(CI)
}

#Orthonormal Series Variance Estimator
SeriesHACFast = function(Preg, Postg, eta, h) {
  library(comprehenr)
  library(stats)
  
  K = nrow(Preg) + 1
  T0 = ncol(Preg)
  T1 = ncol(Postg)
  BarPreg = rowMeans(Preg)
  BarPostg = mean(Postg)
  for (k in nrow(Preg)) {
    Preg[k,] = Preg[k,] - BarPreg[k]
  }
  Postg = Postg -BarPostg

  freq = matrix(0, nrow = K, ncol = h)
  
  for (j in 1:(h)) {
    if (j%%2 == 0) {
      for (t in 1:T0) {
        freq[1:(K-1),j] = freq[1:(K-1),j] + sqrt(2)*sinpi(2*j*(t/T0))*Preg[,t]
      }
      
      freq[1:(K-1),j] = freq[1:(K-1),j]/T0
      
      for (t in 1:T1) {
        freq[K,j] = freq[K,j] + sqrt(2)*sinpi(2*j*(t/T1))*Postg[,t]
      }
      
      freq[K,j] = freq[K,j]/T1
      
    } else {
      
      for (t in 1:T0) {
        freq[1:(K-1),j] = freq[1:(K-1),j] + sqrt(2)*cospi(2*j*(t/T0))*Preg[,t]
      }
      
      freq[1:(K-1),j] = freq[1:(K-1),j]/T0
      
      for (t in 1:T1) {
        freq[K,j] = freq[K,j] + sqrt(2)*cospi(2*j*(t/T1))*Postg[,t]
      }
      freq[K,j] = freq[K,j]/T1
    }
  }
  Vg =  freq %*% t(freq)/h

  V = min(c(T0,T1))*t(eta) %*% Vg %*% eta/eta[K]^2

  return(V)
}


#Find CPE-optimal smoothing parameter using method of Sun (2013)
CPEOptimalh = function(Preg, Postg, p = 1, sig=.05) {
  library(pracma)
  #p: number of parameters of interest
  #sig: significance level
  
  delta2 = qchisq(.75, df = p)
  cv <- qchisq(1 - sig, df = 1)
  cva <- cv       
  tao <- 1.15     
  v <- Preg
  T_val <- ncol(v)
  d <- nrow(v)
  
  # VAR(1) plug in
  dep <- v[, 2:T_val]               
  dep <- t(dep)                     
  indep <- v[, 1:(T_val - 1)]       
  indep <- t(indep)                
    Iindep <- corpcor::pseudoinverse(t(indep) %*% indep) %*% t(indep)
  A <- (t(dep) %*% indep) %*% corpcor::pseudoinverse(t(indep) %*% indep) 
  res <- t(dep) - A %*% t(indep)
  nr <- 1 #number of regressors
  VA <- (res %*% t(res)) / (T_val - nr)
  IA <- corpcor::pseudoinverse(diag(d) - A)
  
  # plug-in estimate of the LRV
  omega0 <- IA %*% VA %*% t(IA)
  
  
  #omegaq1 calculation
  temp <- corpcor::pseudoinverse(diag(nrow(A)) - A)
  part1 <- A %*% VA
  part2 <- A %*% A %*% VA %*% t(A)
  part3 <- A %*% A %*% VA
  part4 <- -6 * A %*% VA %*% t(A)
  part5 <- VA %*% t(A) %*% A
  part6 <- A %*% VA %*% t(A) %*% t(A)
  part7 <- VA %*% t(A)
  
  omegaq1 <- temp %*% temp %*% temp %*% 
    (part1 + part2 + part3 + part4 + part5 + part6 + part7) %*%
    t(temp) %*% t(temp) %*% t(temp)
  
  #plug-in estimate of B
  MB <- - (pi^2) / 6 * omegaq1
  MBbar <- sum(diag(MB %*% corpcor::pseudoinverse(omega0))) / d
  
  # Define density functions
  nchi2den <- function(df, nonc, x) {
    dchisq(x, df = df, ncp = nonc)
  }
  chi2den <- function(df, x) {
    dchisq(x, df = df)
  }
  
  if (MBbar > 0) {
    a1 <- 4 * nchi2den(d, delta2, cva) * abs(MBbar)
    a2 <- delta2 * nchi2den(d + 2, delta2, cva) * 1
    a <- a2 / a1
    Ktemp <- (a^(1/3)) * (T_val^(2/3))
    
  } else if (MBbar <= 0) {
    a1 <- chi2den(d, cva) * cva * abs(MBbar)
    a2 <- (tao - 1) * (1 - sig)
    a <- a2 / a1
    Ktemp <- (a^(1/2)) * T_val
  }
  
  if (Ktemp <= nr + 4) {
    Kwtemp <- nr + 4
  } else if (Ktemp > nr + 4 & Ktemp <= T_val) {
    Kwtemp <- Ktemp
  } else if (Ktemp > T_val) {
    Kwtemp <- T_val
  }
  
  maxkstar <- max(floor(Kwtemp / 2), 1)
  kstar <- max(maxkstar)
  k <- 2 * kstar
  return(k)
}  
