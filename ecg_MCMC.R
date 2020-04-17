proposal_sd <- 0.5 #assume same for alpha_curr, beta, delta_curr, DC for now (the paper did not specify proposal standard deviation ie stepsize)

# this proposal pdf function is in Jordan's lab code but I think we dont need to calculate this in this case cause we assume symmetric proposal 
proposal_lpdf <- function(theta_curr){
  lpdf <- dnorm(theta_curr, proposal_sd)
  return(lpdf)
}

# prior pdf for all model parameters (used gamma(0.01,0.01) for uninformative prior)
prior_lpdf <- function(theta_star, theta_curr){
  lpdf <- dgamma(theta_star, 0.01, 0.01, log = T) - dgamma(theta_curr, 0.01, 0.01, log = T)
  return(lpdf)
}


# likelihood pdf for all model parameters
likelihood_lpdf <- function(y, mean_star, mean_curr, tao_curr){
  lpdf <- 0
  K <- length(y)
  for(k in 1:K){
    lpdf <- lpdf + sum( dnorm(y[k], mean = mean_star[k], sd = tao_curr, log = T) ) - sum ( dnorm(y[k], mean = mean_curr[k], sd = tao_curr, log = T) )
  }
  return(lpdf)
}

# this function calculates the new miu based on the updated theta, theta_star
gen_miu_star <- function(t, alpha_curr, beta_curr, delta_curr) {
  miu_star <- rep(0, 1800)
  for (i in 1:length(t)){
    if (0 < t[i] && t[i] < delta_curr[1]) {
      miu_star[i] <- alpha_curr[1] * cos(pi * t[i] / delta_curr[1]) + beta_curr
    }
    
    
    if (delta_curr[1] < t[i] && t[i] < delta_curr[2]) {
      miu_star[i] <- alpha_curr[2] * (1 - cos(pi * t[i] * (t[i] - delta_curr[1]) / (delta_curr[2]-delta_curr[1]))) + beta_curr - alpha_curr[1]
    }
    
    
    if (delta_curr[2] < t[i] && t[i] < delta_curr[3]) {
      miu_star[i] <- alpha_curr[3] * (1 - cos(pi * (t[i] - delta_curr[2]) / (delta_curr[3]-delta_curr[2]))) + beta_curr - alpha_curr[1] + 2 * alpha_curr[2]
    }
    
    
    if (delta_curr[3] < t[i] && t[i] < delta_curr[4]) {
      miu_star[i] <- 0.5 * (beta_curr - alpha_curr[1] + 2 * alpha_curr[2] + 2 * alpha_curr[3]) * (1 + cos(pi * (t[i] - delta_curr[3]) / (delta_curr[4]-delta_curr[3])))
    }
    
    if (delta_curr[4] < t[i] && t[i] < delta_curr[5]) {
      miu_star[i] <- 0
    }
    
    if (delta_curr[5] < t[i] && t[i]< delta_curr[6]) {
      miu_star[i] <- alpha_curr[4] * (1 - cos(pi * (t[i] - delta_curr[5]) / (delta_curr[6]-delta_curr[5]))) 
    }
    
    
    if (delta_curr[6] < t[i] && t[i] < delta_curr[7]) {
      miu_star[i] <- alpha_curr[4] * (1 + cos(pi * (t[i] - delta_curr[6]) / (delta_curr[7]-delta_curr[6]))) 
    }
    
    
    if (delta_curr[7] < t[i] && t[i] < delta_curr[8]) {
      miu_star[i] <- alpha_curr[5] * (cos(pi * (t[i] - delta_curr[7]) / (delta_curr[8] - delta_curr[7])) - 1) 
    }
    
    
    if (delta_curr[8] < t[i] && t[i] < delta_curr[9]) {
      miu_star[i] <- alpha_curr[6] * (1 - cos(pi * (t[i] - delta_curr[8]) / (delta_curr[9] - delta_curr[8]))) - 2 * alpha_curr[5]
    }
    
    
    if (delta_curr[9] < t[i] && t[i] < delta_curr[10]) {
      miu_star[i] <- alpha_curr[7] * (cos(pi * (t[i] - delta_curr[9]) / (delta_curr[10] - delta_curr[9])) - 1)  + 2 * alpha_curr[6] - 2 * alpha_curr[5]
    }
    
    
    if (delta_curr[10] < t[i] && t[i] < delta_curr[11]) {
      miu_star[i] <- alpha_curr[8] * (1 - cos(pi * (t[i] - delta_curr[10]) / (delta_curr[11] - delta_curr[10]))) - 2 * alpha_curr[7] + 2 * alpha_curr[6] - 2 * alpha_curr[5]
    }
    
    if (delta_curr[11] < t[i] && t[i] < delta_curr[12]) {
      miu_star[i] <- alpha_curr[9] * (1 - cos(pi * (t[i] - delta_curr[11]) / (delta_curr[12] - delta_curr[11]))) + 2 * alpha_curr[8] - alpha_curr[7] + 2 * alpha_curr[6] - 2 * alpha_curr[5]
    }
    
    
    if (delta_curr[12] < t[i] && t[i] < delta_curr[13]) {
      miu_star[i] <- (alpha_curr[9] + alpha_curr[8] - alpha_curr[7] + alpha_curr[6] - alpha_curr[5]) * (1 + cos(pi * (t[i] - delta_curr[12]) / (delta_curr[13] - delta_curr[12])))
    }
    
    
    if (delta_curr[13] < t[i] && t[i] < delta_curr[14]) {
      miu_star[i] <- (alpha_curr[9] + alpha_curr[8] - alpha_curr[7] + alpha_curr[6] - alpha_curr[5]) * (1 + cos(pi * (t[i] - delta_curr[12]) / (delta_curr[13] - delta_curr[12])))
    } 
    
    
    if (delta_curr[14] < t[i] && t[i] < delta_curr[15]) {
      miu_star[i] <- 0
    } 
    
    
    if (delta_curr[14] < t[i] && t[i] < delta_curr[15]) {
      miu_star[i] <- alpha_curr[10] * (1 - cos(pi * (t[i] - delta_curr[14]) / (delta_curr[15] - delta_curr[14])))
    } 
    
    
    if (delta_curr[15] < t[i] && t[i] < delta_curr[16]) {
      miu_star[i] <- alpha_curr[10] * (1 + cos(pi * (t[i] - delta_curr[15]) / (delta_curr[16] - delta_curr[15])))
    } 
    
    
    if (delta_curr[16] < t[i] && t[i] < delta_curr[17]) {
      miu_star[i] <- alpha_curr[11] * (cos(pi * (t[i] - delta_curr[16]) / (delta_curr[17] - delta_curr[16])) - 1)
    } 
    
    
    if (delta_curr[17] < t[i] && t[i] < end_time) {
      miu_star[i] <- alpha_curr[12] * (1 - cos(pi * (t[i] - delta_curr[17]) / (end_time - delta_curr[17]))) - 2 * alpha_curr[11]
    } 
  }
  
  return (miu_star)
}

#function for generating samples for beta

beta_one_samp <- function(beta_curr, BETA, s) { 
  #generate candidate value
  beta_star <- rnorm(1, beta_curr, proposal_sd)

  #
  miu_star <- gen_miu_star(t, alpha_curr = alpha_curr, beta_curr = beta_star, delta_curr = delta_curr)

  #
  mean_star <- miu_star + dc_curr 

  #
  r <- exp((likelihood_lpdf(y, mean_star, mean_curr, tao_curr) + prior_lpdf(beta_star, beta_curr)))
  #print(r)  
  if(runif(1) < r){
    if(s > burn){
      beta_curr <- beta_star
    }
  }
  BETA <- c(BETA, beta_curr)
  return (list(beta_curr, BETA))
}


#function for generating samples for alphai

alpha_one_samp <- function(i, alpha_i_curr, ALPHA_i, s) { 
  #generate candidate value
  alpha_i_star <- rnorm(1, alpha_i_curr, proposal_sd)
  
  #
  alpha_proposed <- alpha_curr
  alpha_proposed[[k]] <- alpha_i_star #alpha_proposed is pointing to a different object than alpha_curr; no change to the list alpha_curr in the global is made
  miu_star <- gen_miu_star(t, alpha_curr = alpha_proposed, beta_curr = beta_curr, delta_curr = delta_curr)
  
  #
  mean_star <- miu_star + dc_curr 
  
  #
  r <- exp((likelihood_lpdf(y, mean_star, mean_curr, tao_curr) + prior_lpdf(alpha_i_star, alpha_i_curr)))
  
  if(runif(1) < r){
    if(s > burn){
      alpha_i_curr <- alpha_i_star
    }
  }
  ALPHA_i <- c(ALPHA_i, alpha_i_curr)
  return (list(alpha_i_curr, ALPHA_i))
}

#function for generating samples for deltai

delta_one_samp <- function(delta_i_curr, delta_i_before, delta_i_next, DELTA_i, s) { 
  #generate candidate value
  delta_i_star <- rnorm(1, delta_i_curr, proposal_sd)
  
  #delta i proposed must be between delta i-1 and next delta i+1, if not just reject it
  if (delta_i_star > delta_i_before && delta_i_star < delta_i_next) { 
    #
    delta_proposed <- delta_curr
    delta_proposed[[k]] <- delta_i_star #delta_proposed is pointing to a different object than delta_curr; no change to the list delta_curr in the global is made
    miu_star <- gen_miu_star(t, alpha_curr = alpha_curr, beta_curr = beta_curr, delta_curr = delta_proposed)
    
    #
    mean_star <- miu_star + dc_curr 
    
    #
    r <- exp((likelihood_lpdf(y, mean_star, mean_curr, tao_curr) + prior_lpdf(delta_i_star, delta_i_curr)))
    
    if(runif(1) < r){
      if(s > burn){
        delta_i_curr <- delta_i_star
      }
    }
  } 
  
  DELTA_i <- c(DELTA_i, delta_i_curr)
  return (list(delta_i_curr, DELTA_i))
}


#function for generating samples for tao

tao_one_samp <- function(tao_curr, BETA, s) { 
  #generate candidate value
  tao_star <- rgamma(1, 0.01, 0.01)
  
  #
  miu_star <- gen_miu_star(t, alpha_curr = alpha_curr, beta_curr = beta_curr, delta_curr = delta_curr)
  
  #
  mean_star <- miu_star + dc_curr 
  
  # note that tao in likelihood_lpdf is now tao_star
  r <- exp((likelihood_lpdf(y, mean_star, mean_curr, tao_star) + prior_lpdf(tao_star, tao_curr)))
  #print(r)  
  if(runif(1) < r){
    if(s > burn){
      tao_curr <- tao_star
    }
  }
  TAO <- c(TAO, tao_curr)
  return (list(tao_curr, TAO))
}

#function for generating samples for DC

dc_one_samp <- function(dc_curr, DC, s) { 
  #generate candidate value
  dc_star <- rnorm(1, dc_curr, proposal_sd)
  
  #
  miu_star <- gen_miu_star(t, alpha_curr = alpha_curr, beta_curr = beta_curr, delta_curr = delta_curr)
  
  # note that dc_star is now used instead of dc_curr
  mean_star <- miu_star + dc_star 
  
  #
  r <- exp((likelihood_lpdf(y, mean_star, mean_curr, tao_curr) + prior_lpdf(dc_star, dc_curr)))
  #print(r)  
  if(runif(1) < r){
    if(s > burn){
      dc_curr <- dc_star
    }
  }
  DC <- c(DC, dc_curr)
  return (list(dc_curr, DC))
}


#Initiate 

hr_data_1 <- read.table("~/projects/STA360-Final_Project/hr_data_1.11839", quote="\"", comment.char="")
S <- 10000
burn <- 5000
y <- hr_data_1$V1
t <- seq(0, 900, by = 0.5)
beta_curr <- 1
alpha_curr <- rep(1, 12)
delta_curr<- seq(20, 860, length.out = 17) 
end_time <- 900
tao_curr <- 0.0001
mean_curr <- rep(mean(y), length(y)) #mean_curr is a vector cause mean_star of y is a vector so I think r only lets me run it if it is a vector
dc_curr <- 5
miu_curr <- mean_curr - dc_curr

#storage

DC <- NULL
ALPHA <- rep(list(NULL), 12)
BETA <- NULL
DELTA <- rep(list(NULL), 17)
TAO <- NULL

for(s in 1:S) {
  
  #1 draw of Beta
  
  beta_lst <- beta_one_samp(beta_curr, BETA, s)
  beta_curr <- beta_lst[[1]]
  BETA <- beta_lst[[2]]
  
  
  #1 draw of Alphas (1-12)
  
  for (k in 1:12) {
    alpha_lst <- alpha_one_samp(k, alpha_curr[k], ALPHA[[k]], s)
    alpha_curr[k] <- alpha_lst[[1]]
    ALPHA[[k]] <- alpha_lst[[2]]
  }
  
  #1 draw of Deltas (1-17)
  for (k in 1:17) {
    delta_lst <- delta_one_samp(k, delta_curr[k], DELTA[[k]], s)
    delta_curr[k] <- delta_lst[[1]]
    DELTA[[k]] <- delta_lst[[2]]
  }
  
  #1 draw of Tao 
  tao_lst <- tao_one_samp(tao_curr, TAO, s)
  tao_curr <- tao_lst[[1]]
  TAO <- tao_lst[[2]]
  
  
  #1 draw of DC 
  dc_lst <- dc_one_samp(dc_curr, DC, s)
  dc_curr <- dc_lst[[1]]
  DC <- dc_lst[[2]]
  
}