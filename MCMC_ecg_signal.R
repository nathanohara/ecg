proposal_sd <- 0.5 #assume same for alpha, beta, delta, DC for now
prior_sd <- 0.5 #assume same for alpha, beta, delta, DC for now

proposal_lpdf <- function(theta_curr){
  lpdf <- dnorm(theta_curr, proposal_sd)
  return(lpdf)
}

#
prior_lpdf_ <- function(theta_star, theta_curr){
  lpdf <- lpdf + dnorm(theta_star, sd = prior_sd, log = T) - dnorm(theta_curr, sd = prior_sd, log = T)
  return(lpdf)
}

#
likelihood_lpdf <- function(y, mean_star, mean_curr, tao){
  lpdf <- 0
  K <- length(y)
  for(k in 1:K){
    lpdf <- lpdf + sum( dnorm(y[k], mean = mean_star[k], sd = tao, log = T) ) - sum ( dnorm(y[k], mean = mean_curr[k], sd = tao, log = T) )
  }
  return(lpdf)
}

#
gen_miu_star <- function(delta_curr) {
  if (0 < delta_curr & delta_curr < delta[1]) {
    miu_star <- alpha[1] * cos(pi * t / delta[1]) + beta_param
  }
  
  
  if (delta[1] < delta_curr & delta_curr < delta[2]) {
    miu_star <- alpha[2] * (1 - cos(pi * t * (t - delta[1]) / (delta[2]-delta[1]))) + beta_param - alpha[1]
  }
  
  
  if (delta[2] < delta_curr & delta_curr < delta[3]) {
    miu_star <- alpha[3] * (1 - cos(pi * (t - delta[2]) / (delta[3]-delta[2]))) + beta_param - alpha[1] + 2 * alpha[2]
  }
  
  
  if (delta[3] < delta_curr & delta_curr < delta[4]) {
    miu_star <- 0.5 * (beta_param - alpha[1] + 2 * alpha[2] + 2 * alpha[3]) * (1 + cos(pi * (t - delta[3]) / (delta[4]-delta[3])))
  }
  
  if (delta[4] < delta_curr & delta_curr < delta[5]) {
    miu_star <- 0
  }
  
  if (delta[5] < delta_curr & delta_curr< delta[6]) {
    miu_star <- alpha[4] * (1 - cos(pi * (t - delta[5]) / (delta[6]-delta[5]))) 
  }
  
  
  if (delta[6] < delta_curr & delta_curr < delta[7]) {
    miu_star <- alpha[4] * (1 + cos(pi * (t - delta[6]) / (delta[7]-delta[6]))) 
  }
  
  
  if (delta[7] < delta_curr & delta_curr < delta[8]) {
    miu_star <- alpha[5] * (cos(pi * (t - delta[7]) / (delta[8] - delta[7])) - 1) 
  }
  
  
  if (delta[8] < delta_curr & delta_curr < delta[9]) {
    miu_star <- alpha[6] * (1 - cos(pi * (t - delta[8]) / (delta[9] - delta[8]))) - 2 * alpha[5]
  }
  
  
  if (delta[9] < delta_curr & delta_curr < delta[10]) {
    miu_star <- alpha[7] * (cos(pi * (t - delta[9]) / (delta[10] - delta[9])) - 1)  + 2 * alpha[6] - 2 * alpha[5]
  }
  
  
  if (delta[10] < delta_curr & delta_curr < delta[11]) {
    miu_star <- alpha[8] * (1 - cos(pi * (t - delta[10]) / (delta[11] - delta[10]))) - 2 * alpha[7] + 2 * alpha[6] - 2 * alpha[5]
  }
  
  if (delta[11] < delta_curr & delta_curr < delta[12]) {
    miu_star <- alpha[9] * (1 - cos(pi * (t - delta[11]) / (delta[12] - delta[11]))) + 2 * alpha[8] - alpha[7] + 2 * alpha[6] - 2 * alpha[5]
  }
  
  
  if (delta[12] < delta_curr & delta_curr < delta[13]) {
    miu_star <- (alpha[9] + alpha[8] - alpha[7] + alpha[6] - alpha[5]) * (1 + cos(pi * (t - delta[12]) / (delta[13] - delta[12])))
  }
  
  
  if (delta[13] < delta_curr & delta_curr < delta[14]) {
    miu_star <- (alpha[9] + alpha[8] - alpha[7] + alpha[6] - alpha[5]) * (1 + cos(pi * (t - delta[12]) / (delta[13] - delta[12])))
  } 
  
  
  if (delta[14] < delta_curr & delta_curr < delta[15]) {
    miu_star <- 0
  } 
  
  
  if (delta[14] < delta_curr & delta_curr < delta[15]) {
    miu_star <- alpha[10] * (1 - cos(pi * (t - delta[14]) / (delta[15] - delta[14])))
  } 
  
  
  if (delta[15] < delta_curr & delta_curr < delta[16]) {
    miu_star <- alpha[10] * (1 + cos(pi * (t - delta[15]) / (delta[16] - delta[15])))
  } 
  
  
  if (delta[16] < delta_curr & delta_curr < delta[17]) {
    miu_star <- alpha[11] * (cos(pi * (t - delta[16]) / (delta[17] - delta[16])) - 1)
  } 
  
  
  if (delta[17] < delta_curr & delta_curr < T) {
    miu_star <- alpha[12] * (1 - cos(pi * (t - delta[17]) / (T - delta[17]))) - 2 * alpha[11]
  } 
  return (miu_star)
}

#function for generating sampels for 1 model parameter

one_param <- function(theta_curr, THETA, s, burn, proposal_sd, mean_curr, dc_curr, tao) { #when dc is the parameter of interest, theta_curr = dc_curr
  #generate candidate value
  theta_star <- rnorm(1, theta_curr, proposal_sd)
    
  #
  miu_star <- gen_miu_star(theta_star)
    
  #
  mean_star <- miu_star + dc_curr 
    
  #
  r <- exp((likelihood_lpdf(y, mean_star, mean_curr, tao) + prior_lpdf(theta_star, theta_curr)))
    
  if(runif(1) < r){
    if(s > burn){
      theta_curr <- theta_star
    }
  }
  THETA <- c(THETA, theta_curr)
}



S <- 10000
burn <- 5000
y <- hr_data_1$V1
t <- seq(0, 900, by = 0.5)
beta_param <- 1
alpha <- rep(1, 12)
delta <- seq(20, 860, length.out = 17) 
T <- 900
mean_curr <- mean(y)
dc_curr <- 5
miu_curr <- mean_curr - dc_curr

#storage

MIU <- NULL
DC <- NULL
ALPHA <- NULL
BETA <- NULL
DETA <- NULL
