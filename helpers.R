
#
# Weibull ####

rate2scale <- function(rate, shape){
  rate^(-1/shape)
}

rweibull2 <- function(n, shape, rate){
  scale   <- rate2scale(rate, shape)
  rweibull(n, shape = shape, scale = scale)
}

dweibull2 <- function(x, shape, rate, log){
  scale   <- rate2scale(rate, shape)
  dweibull(x, shape = shape, scale = scale, log = log)
}

pweibull2 <- function(x, shape, rate, log, lower.tail){
  scale   <- rate2scale(rate, shape)
  pweibull(x, shape = shape, scale = scale, log.p = log, lower.tail = lower.tail)
}

pweibull2_ab <- function(a, b, shape, rate, log, lower.tail){
  pweibull2(b, shape, rate, log = F, lower.tail = F) - 
  pweibull2(a, shape, rate, log = F, lower.tail = F)
}


#
# Process ####

draw_sojourn_H <- function(theta){
  rweibull2(1, shape = theta$shape_H, rate = theta$rate_H)
}

draw_sojourn_P <- function(theta){
  rweibull2(1, shape = theta$shape_P, rate = theta$rate_P)
}

compute_compartment <- function(t, tau_HP, tau_PC){
  case_when(
    t < tau_HP ~ "H",
    TRUE       ~ "P"
  )
}

screen_result <- function(compartment, theta) {
  if(compartment == "H")  return(FALSE)  
  if(compartment == "P")  return(rbernoulli(1, p = theta$beta))
}

#
# Gibbs ####

gibbs_rate_H <- function(tau_HP, theta, prior){
  
  a_n <- prior$shape_H + length(tau_HP)
  
  sojourn_H   <- tau_HP - 40
  sum_sojourn <- sum(sojourn_H^theta$shape_H) 
  b_n         <- prior$rate_H + sum_sojourn
  
  rgamma(1, shape = a_n, rate = b_n)
}

gibbs_rate_P <- function(d_obs_censor, tau_HP, indolent, theta, prior){
  
  a_n <- prior$shape_P + sum(d_obs_censor$censor_type == "clinical")
  
  sojourn_P   <- (d_obs_censor$censor_time - tau_HP)[!indolent]
  sum_sojourn <- sum(sojourn_P^theta$shape_P) 
  b_n         <- prior$rate_P + sum_sojourn
  
  rgamma(1, shape = a_n, rate = b_n)
}

gibbs_psy <- function(indolent, prior){
  
  a_n <- prior$a_psy + sum(indolent  ) 
  b_n <- prior$b_psy + sum(1-indolent) 
  
  rbeta(1, a_n, b_n)
}

gibbs_beta <- function(d_process, d_obs_screen, prior){
  
  screen_preclinical <- d_obs_screen %>%
    left_join(d_process, by = "person_id") %>%
    mutate(compartment = compute_compartment(age_screen, tau_HP, tau_PC)) %>%
    filter(compartment == "P") %>%
    pull(screen_detected)
  
  a_0 <- prior$a_beta
  x   <- sum(screen_preclinical)
  a_n <- a_0 + x
  
  b_0 <- prior$b_beta
  n   <- length(screen_preclinical)
  b_n <- b_0 + (n-x)
  
  rbeta(1, a_n, b_n)
}

#
# Likelihood ####

dlog_sojourn_H <- function(sojourn_H, theta){
  dweibull2(
    sojourn_H, 
    shape = theta$shape_H, rate = theta$rate_H,
    log = TRUE
  )
}

dlog_sojourn_P <- function(sojourn_P, theta, indolent, censor_type){
  
  if(indolent){
    0
  }else if(censor_type == "clinical"){
    dweibull2(
      sojourn_P, 
      shape = theta$shape_P, rate = theta$rate_P, 
      log = TRUE
    )
  } else if(censor_type == "screen"){
    pweibull2(
      sojourn_P, 
      shape = theta$shape_P, rate = theta$rate_P, 
      log = TRUE, lower.tail = FALSE
    )
  }
  
}

dlog_psy <- function(indolent, theta){
  dbinom(
    indolent, 1, 
    prob = theta$psy, log = TRUE
  )
}

dlog_beta <- function(d_obs_screen_i, theta, tau_HP, censor_time){
  
  # Extract outcome of screens during preclinical phase
  screen_preclinical <- d_obs_screen_i %>%
    mutate(compartment = compute_compartment(age_screen, tau_HP, censor_time)) %>%
    filter(compartment == "P") %>%
    pull(screen_detected)
  
  n_screen          <- length(screen_preclinical)
  n_screen_positive <- sum(screen_preclinical)
  
  # product of Bernoulli RVs
  n_screen_positive*log(theta$beta) + (n_screen - n_screen_positive)*log(1-theta$beta)
  
}

dlog_likelihood_i <- function(tau_HP, indolent, d_obs_screen_i, d_obs_censor_i, theta){
  
  censor_time <- d_obs_censor_i$censor_time
  
  dlog_H <- dlog_sojourn_H(tau_HP - 40         , theta)
  dlog_P <- dlog_sojourn_P(censor_time - tau_HP, theta, indolent, d_obs_censor_i$censor_type)
  dlog_p <- dlog_psy      (indolent, theta)
  dlog_b <- dlog_beta     (d_obs_screen_i, theta, tau_HP, censor_time)
  
  dlog_likelihood <- dlog_H + dlog_P + dlog_p + dlog_b
  return(dlog_likelihood)
  
}


#
# MH ####
rprop_tau_HP <- function(d_obs_screen_i, d_obs_censor_i, theta){
  
  
  # construct interval (t_i, t_{i+1}]
  age_screen <- d_obs_screen_i$age_screen # age at screens
  if(d_obs_censor_i$censor_type == "clinical"){ # if clinical, add tau_PC as endpoint of last interval
    age_screen <- c(age_screen, d_obs_censor_i$censor_time)
  }
  dt         <- diff(age_screen) # length of intervals
  
  # construct probability vector p
  K          <- length(dt)
  prob       <- numeric(K)
  
  for(k in 1:K){
    prob[k] <- if(d_obs_censor_i$censor_type == "screen"  ){
          theta$beta * (1-theta$beta)^(K-k) * dt[k]
        } else if(d_obs_censor_i$censor_type == "clinical"){
          (1-theta$beta)^(K-k) * dt[k]
        } # end-if
  } # end-for
  
  prob <- prob/sum(prob) # normalize probabilities
  
  # sample tau_HP
  k_new      <- sample.int(K, 1, prob = prob) # sample interval from {1, ..., K}
  tau_HP_new <- runif(1, age_screen[k_new], age_screen[k_new+1]) # sample tau_HP in the chosen interval
  
  return(tau_HP_new)

}

dlog_prop_tau_HP <- function(tau_HP, d_obs_screen_i, d_obs_censor_i, theta){
  
  # same code as in rprop_tau_HP
  {
    # construct interval (t_i, t_{i+1}]
    age_screen <- d_obs_screen_i$age_screen # age at screens
    if(d_obs_censor_i$censor_type == "clinical"){ # if clinical, add tau_PC as endpoint of last interval
      age_screen <- c(age_screen, d_obs_censor_i$censor_time)
    }
    dt         <- diff(age_screen) # length of intervals
    
    # construct probability vector p
    K          <- length(dt)
    prob       <- numeric(K)
    
    for(k in 1:K){
      prob[k] <- if(d_obs_censor_i$censor_type == "screen"  ){
        theta$beta * (1-theta$beta)^(K-k) * dt[k]
      } else if(d_obs_censor_i$censor_type == "clinical"){
        (1-theta$beta)^(K-k) * dt[k]
      } # end-if
    } # end-for
    
    prob <- prob/sum(prob) # normalize probabilities
    
  } # end-same code as in rprop_tau_HP
  
  # find k_new
  k_new <- sum(age_screen < tau_HP) 
  
  # log density
  log(prob[k_new]) - log(age_screen[k_new+1] - age_screen[k_new]) # log(categorical(p) * uniform(dt_k))
  
}

MH_tau_PH <- function(i, d_obs_screen, d_obs_censor, indolent, tau_HP, theta){
  
  # data for person i
  d_obs_screen_i <- d_obs_screen %>% filter(person_id == i)
  d_obs_censor_i <- d_obs_censor %>% filter(person_id == i)
  indolent_i     <- indolent[i] 
  
  # propose new tau_HP
  tau_HP_cur <- tau_HP[i] 
  tau_HP_new <- rprop_tau_HP(d_obs_screen_i, d_obs_censor_i, theta)
  
  # M-H acceptance ratio
  dlog_prop_cur  <- dlog_prop_tau_HP (tau_HP_cur, d_obs_screen_i, d_obs_censor_i, theta)
  dlog_prop_new  <- dlog_prop_tau_HP (tau_HP_new, d_obs_screen_i, d_obs_censor_i, theta)
  dlog_lik_cur   <- dlog_likelihood_i(tau_HP_cur, indolent_i, d_obs_screen_i, d_obs_censor_i, theta)
  dlog_lik_new   <- dlog_likelihood_i(tau_HP_new, indolent_i, d_obs_screen_i, d_obs_censor_i, theta)

  MH_logratio <- dlog_lik_new - dlog_lik_cur + dlog_prop_cur - dlog_prop_new
  
  if(runif(1) < exp(MH_logratio)){ # accept new tau_HP
    tau_HP_cur <- tau_HP_new
  }
  
  return(tau_HP_cur)
  
}

#
# MCMC ####
MCMC <- function(){
  
  { # runtime
  m <- 1e2
  RATE_H <- RATE_P <- BETA <- PSY <- numeric(m)
  indolent <- d_process$indolent
  tau_HP_tmp <- tau_HP   <- d_process$tau_HP
  theta_tmp <- theta
  
  tic()
  for(n in 1 : m){ print(n)
    
    # update parameters
    rate_H    <- gibbs_rate_H(tau_HP, theta_tmp, prior)
    rate_P    <- gibbs_rate_P(d_obs_censor, tau_HP, indolent, theta_tmp, prior)
    psy       <- gibbs_psy   (indolent , prior)
    beta      <- gibbs_beta  (d_process, d_obs_screen, prior)
    theta_tmp <- list(rate_H=rate_H, rate_P=rate_P, psy=psy, beta=beta,
                      shape_H=theta_tmp$shape_H, shape_P=theta_tmp$shape_P)
    
    # update latent data
    for(i in 1 : nrow(d_process)){ # this for-loop should be done in parallel
      tau_HP[i] <- MH_tau_PH(i, d_obs_screen, d_obs_censor, indolent, tau_HP, theta_tmp)
    }
    
    print(mean(tau_HP_tmp != tau_HP)) # acceptance rate
    tau_HP_tmp <- tau_HP
    
    RATE_H[n] <- rate_H
    RATE_P[n] <- rate_P
    BETA  [n] <- beta
    PSY   [n] <- psy
    
  }
  runtime <- toc()
  } # end runtime
  
  # output
  THETA <- list(RATE_H=RATE_H, RATE_P=RATE_P, BETA=BETA, PSY=PSY)
  out <- list(THETA)
  return(out)
  
}
