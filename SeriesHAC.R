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
  
  for (j in 0:(h-1)) {
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

SeriesHAC = function(Preg, Postg, eta, h, alpha) {
  library(comprehenr)
  
  K = nrow(Preg) + 1
  T0 = ncol(Preg)
  T1 = ncol(Postg)
  BarPreg = rowMeans(Preg)
  BarPostg = mean(Postg)
  #Tmax = max(c(T0,T1))
  #g = matrix(0, nrow = K, ncol = Tmax)
  #for (k in 1:(K-1)) {
  #  g[k,1:T0] = Preg[k,] - BarPreg[k]
  #}
  #g[K+1,1:T1] = Postg - BarPostg 
  
  #Vg = sum(to_vec(for (s in 1:Tmax) (sum(to_vec(for (i in 1:Tmax) Q(i/Tmax,s/Tmax,h) as.matrix(g[,s]) %*% t(as.matrix(g[,i])) )))))
  
  Vg = matrix(NA, nrow = K, ncol = K)
  for (g in 1:K) {
    for (l in 1:K) {
      if (g < K & l < K) {
        Vg[g,l] = sum(to_vec(for (s in 1:T0) (sum(to_vec(for (i in 1:T0) Q(i/T0,s/T0,h)*(Preg[g,i] - BarPreg[g])*(Preg[l,s] - BarPreg[l]))))))/T0^2
      }
      if (g == K & l < K) {
        Vg[g,l] = sum(to_vec(for (s in 1:T0) (sum(to_vec(for (i in 1:T1) Q(i/T1,s/T0,h)*(Postg[1,i] - BarPostg)*(Preg[l,s] - BarPreg[l]))))))/(T0*T1)
      }
      if (g < K & l == K) {
        Vg[g,l] = sum(to_vec(for (s in 1:T1) (sum(to_vec(for (i in 1:T0) Q(i/T0,s/T1,h)*(Preg[g,i] - BarPreg[g])*(Postg[1,s] - BarPostg))))))/(T0*T1)
      }
      if (g == K & l == K) {
        Vg[g,l] = sum(to_vec(for (s in 1:T1) (sum(to_vec(for (i in 1:T1) Q(i/T1,s/T1,h)*(Postg[1,i] - BarPostg)*(Postg[1,s] - BarPostg))))))/T1^2
      }
    }
  }
  V = min(c(T0,T1))*t(eta) %*% Vg %*% eta/eta[K]^2
  return(V)
}

Optimalh = function(Preg, Postg, p, alpha) {
  library(vars)
  
  #initialVg = SeriesHAC(Preg, Postg, initialh)
  K = nrow(Preg) + 1
  T0 = ncol(Preg)
  T1 = ncol(Postg)
  n = min(c(T0,T1)) 
  
  #Get initial Var(1) estimate
  var1 = VAR(t(Preg), p = 1)
  B = BQ(var1)$B
  Omega0 = cov(t(Preg))
  
  MB=-(pi^2)/6*B
  MBbar=sum(diag((MB %*% inv(Omega0))))/p
  
  
  delta2 = qchisq(.75, df = p)
  d = qchisq(1-alpha, df = p)
  a1 = dchisq(d, df = p+2, ncp = delta2)
  a2 = dchisq(d, df = p, ncp = delta2)
  a3 = dchisq(d, df = p)
  tau = 2
  if (MBbar > 0) {
    Kwtemp = ((delta2*a1)/(4*a2*abs(MBbar)))^(1/3)*n^(2/3)
  } else {
    Kwtemp = (((tau-1)*alpha)/(a3*d*abs(MBbar)))^(1/2)*n
  }
  
  if (Kwtemp < (p+4)) {
    Kwtemp= p+4
  }
  if(Kwtemp > n) {
    Kwtemp = n
  }
  
  Kwtemp = floor(Kwtemp/2)
  h = 2*Kwtemp
  
  return(h)
}