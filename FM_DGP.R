TWFE_DGP <- function(J = NULL, T = NULL, unitFE = NULL, timeFE_par = NA, Sigma = NULL, timeFE_estimated = FALSE) {
  #Initializing Matrix
  data <- as.data.frame(matrix(NA, J*T,4))
  colnames(data) <- c("unit.num","year","Y","name")
  # Create the unit identifiers.
  for (j in 1:J) {
    data[(1 + (j-1)*T):(T + (j-1)*T), 1] <- rep(j, T)
  }
  #Creating unit names
  data[,4] <- as.character(data[,1])
  # Create the time identifiers.
  data[, 2] <- rep(seq(from = 1, to  = T, by = 1), J)
  # Create time fixed effects
  if (!timeFE_estimated) {
      rho <- 0.9
      sd.rho <- sqrt(1-rho^2)
      delta <- matrix(arima.sim(list(ar = rho), n = T, sd=sd.rho))
  } else {
    delta = forecast(timeFE_par, h = T)$mean
  }
  # Create unit fixed effects
  mu <- unitFE
  #Creating untreated and treated outcomes
  for (j in 1:J) {
    data[(1 + (j-1)*T):(T + (j-1)*T),3] <- delta + mu[j] + rnorm(T, mean = 0, sd = Sigma[j])  
  }
  return(data)
}

InterFE_DGP <- function(N = NULL, T = NULL, timeFE_par = NULL, unitFE = NULL,loadings = NULL, Factor_pars = NULL, Sigma = NULL) {
  #Initializing Matrix
  data <- as.data.frame(matrix(NA, N*T,4))
  colnames(data) <- c("unit.num","year","Y","name")
  # Create the unit identifiers.
  for (j in 1:N) {
    data[(1 + (j-1)*T):(T + (j-1)*T), 1] <- rep(j, T)
  }
  #Creating unit names
  data[,4] <- as.character(data[,1])
  # Create the time identifiers.
  data[, 2] <- rep(seq(from = 1, to  = T, by = 1), N)
  #Simulate factors
  F = length(Factor_pars)
  Factors = matrix(NA, T, F)
  for(f in 1:F) {
    Factors[,f] = forecast(Factor_pars[[f]], h = T)$mean
  }
  timeFE = forecast(timeFE_par, h = T)$mean
  #Creating outcomes
  for (j in 1:N) {
    data[(1 + (j-1)*T):(T + (j-1)*T),3] <- timeFE + unitFE[j] + Factors %*% loadings[,j] + rnorm(T, mean = 0, sd = Sigma[j])
  }
  return(data)
}

FM_DGP <- function(N = NULL, Time = NULL, loadings = NULL, Factor_pars = NULL, Sigma = NULL, serial = FALSE) {
  library(VGAM)
  #Initializing Matrix
  data <- as.data.frame(matrix(nrow = N*Time, ncol = 4))
  colnames(data) <- c("Countryno","year","Y","name")
  # Create the unit identifiers.
  for (j in 1:N) {
    data[(1 + (j-1)*Time):(Time + (j-1)*Time), 1] <- rep(j, Time)
  }
  #Creating unit names
  data[,4] <- as.character(data[,1])
  # Create the time identifiers.
  data[, 2] <- rep(seq(from = 1, to  = Time, by = 1), N)
  #Simulate factors
  F = length(Factor_pars)
  Factors = matrix(NA, Time, F)
  for(f in 1:F) {
    Factors[,f] = forecast(Factor_pars[[f]], h = Time)$mean
    #Factors[,f] = predict(Factor_pars[[f]], n.ahead = Time)$pred
    #Factors[,f] = as.numeric(arima.sim(model = Factor_pars[[f]], n = Time))
  }
   
  #Creating outcomes
  for (j in 1:N) {
    #Simulate Shocks
    if (serial) {
      epsilon_j = as.numeric(arima.sim(model = Sigma[[j]], n = Time))
    } else {
      epsilon_j = rnorm(Time, mean = 0, sd = Sigma[j]) #rlaplace(n=T, location=0, scale=Sigma[j]/sqrt(2))
    }
    
    data[(1 + (j-1)*Time):(Time + (j-1)*Time),3] = Factors %*% loadings[,j] + epsilon_j
  }
  return(data)
}

FM_ARMA_DGP <- function(N = NULL, T = NULL, loadings = NULL, Factor_pars = NULL, epsilon_pars = NULL) {
  #Initializing Matrix
  data <- as.data.frame(matrix(NA, N*T,4))
  colnames(data) <- c("unit.num","year","Y","name")
  # Create the unit identifiers.
  for (j in 1:N) {
    data[(1 + (j-1)*T):(T + (j-1)*T), 1] <- rep(j, T)
  }
  #Creating unit names
  data[,4] <- as.character(data[,1])
  # Create the time identifiers.
  data[, 2] <- rep(seq(from = 1, to  = T, by = 1), N)
  #Simulate factors
  F = length(Factor_pars)
  Factors = matrix(NA, T, F)
  for(f in 1:F) {
    Factors[,f] = forecast(Factor_pars[[f]], h = T)$mean
  }
  epsilon = matrix(NA, T, N)
  for(i in 1:N) {
    epsilon[,i] = forecast(epsilon_pars[[i]], h = T)$mean
  }
  #Creating outcomes
  for (j in 1:N) {
    data[(1 + (j-1)*T):(T + (j-1)*T),3] <- Factors %*% loadings[,j] + epsilon[,j]
  }
  return(data)
}
