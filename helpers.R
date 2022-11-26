
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

pweibull2_ab <- function(a, b, shape, rate, log){
  pweibull2(b, shape, rate, log = F, lower.tail = T) - 
  pweibull2(a, shape, rate, log = F, lower.tail = T)
}


#
# Process ####

draw_sojourn_H <- function(theta){
  rweibull2(1, shape = theta$shape_H, rate = theta$rate_H)
}

draw_sojourn_P <- function(theta){
  rweibull2(1, shape = theta$shape_P, rate = theta$rate_P)
}

compute_compartment <- function(t, tau_HP){
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

gibbs_rate_H <- function(tau_HP, theta, prior, d_obs_censor){ # gamma prior and weibull likelihood conjugacy
  
  a_n <- prior$shape_H + sum(is.finite(tau_HP)) # number of observed weibull
  
  sojourn_H   <- pmin(tau_HP, d_obs_censor$censor_time) - 40
  sum_sojourn <- sum(sojourn_H^theta$shape_H) 
  b_n         <- prior$rate_H + sum_sojourn
  
  rgamma(1, shape = a_n, rate = b_n)
}

gibbs_rate_P <- function(tau_HP, indolent, theta, prior, d_obs_censor){ # gamma prior and weibull likelihood conjugacy
  
  a_n <- prior$shape_P + sum(d_obs_censor$censor_type == "clinical") # numbe rof observed weibull
  
  sojourn_P   <- (d_obs_censor$censor_time - tau_HP)[!indolent & is.finite(tau_HP)]
  sum_sojourn <- sum(sojourn_P^theta$shape_P) 
  b_n         <- prior$rate_P + sum_sojourn
  
  rgamma(1, shape = a_n, rate = b_n)
}

gibbs_psy <- function(indolent, prior){ # beta prior and binomial likelihood conjugacy
  
  a_n <- prior$a_psy + sum(indolent  , na.rm = TRUE) 
  b_n <- prior$b_psy + sum(1-indolent, na.rm = TRUE) 
  
  rbeta(1, a_n, b_n)
}

gibbs_beta <- function(tau_HP, d_id, d_obs_screen, prior, n_screen_positive){ # beta prior and binomial likelihood conjugacy
  
  n_screen <- d_id %>% 
    mutate(tau_HP=tau_HP) %>%
    left_join(d_obs_screen, by = "person_id") %>%
    mutate(compartment = compute_compartment(age_screen, tau_HP)) %>%
    filter(compartment == "P") %>%
    nrow()
  
  a_n <- prior$a_beta + n_screen_positive
  b_n <- prior$b_beta + (n_screen-n_screen_positive)
  
  rbeta(1, a_n, b_n)
}

#
# Likelihood ####

dlog_sojourn_H <- function(tau_HP, theta, censor_time){
  
  if(is.infinite(tau_HP)){
    pweibull2(
      censor_time - 40,
      shape = theta$shape_H, rate = theta$rate_H,
      log = TRUE, lower.tail = FALSE
      )
  }else if(is.finite(tau_HP)){
    dweibull2(
      tau_HP - 40, 
      shape = theta$shape_H, rate = theta$rate_H,
      log = TRUE
    )
  }
  
}

dlog_sojourn_P <- function(tau_HP, theta, indolent, censor_type, censor_time){
  
  if(is.infinite(tau_HP)){
    0
  }else if(indolent){
    0
  }else if(censor_type == "clinical"){
    dweibull2(
      censor_time - tau_HP, 
      shape = theta$shape_P, rate = theta$rate_P, 
      log = TRUE
    )
  } else if(censor_type %in% c("screen", "censored")){
    pweibull2(
      censor_time - tau_HP, 
      shape = theta$shape_P, rate = theta$rate_P, 
      log = TRUE, lower.tail = FALSE
    )
  }
  
}

dlog_psy <- function(indolent, theta){
  if(is.na(indolent)){
    0
  } else {
    dbinom(
      indolent, 1, 
      prob = theta$psy, log = TRUE
    )
  }
}

dlog_beta <- function(d_obs_screen, theta, tau_HP, censor_type){
  
  n_screen          <- sum(d_obs_screen$age_screen>tau_HP) # number of screens during pre-clinical  
  n_screen_positive <- censor_type=="screen" # number of positive screens during pre-clinical  
  
  # product of Bernoulli RVs
  n_screen_positive*log(theta$beta) + (n_screen - n_screen_positive)*log(1-theta$beta)
  
}

dlog_likelihood_i <- function(tau_HP_i, indolent_i, d_obs_screen_i, d_obs_censor_i, theta){
  
  censor_time_i <- d_obs_censor_i$censor_time
  censor_type_i <- d_obs_censor_i$censor_type
  
  dlog_H <- dlog_sojourn_H(tau_HP_i, theta, censor_time_i)
  dlog_P <- dlog_sojourn_P(tau_HP_i, theta, indolent_i, censor_type_i, censor_time_i)
  dlog_p <- dlog_psy      (indolent_i, theta)
  dlog_b <- dlog_beta     (d_obs_screen_i, theta, tau_HP_i, censor_type_i)
  
  dlog_likelihood <- dlog_H + dlog_P + dlog_p + dlog_b
  return(dlog_likelihood)
  
}


#
# MH ####
rprop_tau_HP <- function(d_obs_screen_i, d_obs_censor_i, theta){
  
  # same code as in dprop_tau_HP
  {
    censor_type <- d_obs_censor_i$censor_type
    
    # construct interval (t_i, t_{i+1}]
    age_screen <- d_obs_screen_i$age_screen # age at screens
    
    if(censor_type == "clinical"){ # if clinical, add tau_PC as endpoint of last interval
      age_screen <- c(age_screen, d_obs_censor_i$censor_time)
    }else if(censor_type == "censored"){ # if censored, add Inf as endpoint of last interval
      age_screen <- c(age_screen, d_obs_censor_i$censor_time, Inf)
    }
    
    dt         <- diff(age_screen) # length of intervals
    K          <- length(dt)       # number of intervals
    
    # construct probability vector p
    prob <- if(censor_type == "screen"  ){
      
      theta$beta * (1-theta$beta)^((K-1):0) * dt[1:K]
      
    }  else if(censor_type == "clinical"){
      
      (1-theta$beta)^((K-1):0) * dt[1:K]
      
    } else if(censor_type == "censored"){
      
      prob_sojourn_H <- pweibull2_ab(
        age_screen[1:K    ] - 40, 
        age_screen[1:K + 1] - 40, 
        shape = theta$shape_H, rate = theta$rate_H, 
        log = FALSE
      )
      
      prob_screens <- (1-theta$beta)^((K-1):0)
      
      prob_sojourn_H * prob_screens
      
    }
  
    prob <- prob/sum(prob) # normalize probabilities
  }  # end-same code as in dprop_tau_HP
  
  # sample tau_HP
  k_new      <- sample.int(K, 1, prob = prob) # sample interval
  tau_HP_new <- if(censor_type == "censored" & k_new == K){
    Inf # tau_HP > time_censor
  }else{
    runif(1, age_screen[k_new], age_screen[k_new+1]) # sample tau_HP in the chosen interval
  }
  
  return(tau_HP_new)

}

dlog_prop_tau_HP <- function(tau_HP, d_obs_screen_i, d_obs_censor_i, theta){
  
  # same code as in rprop_tau_HP
  {
    censor_type <- d_obs_censor_i$censor_type
    
    # construct interval (t_i, t_{i+1}]
    age_screen <- d_obs_screen_i$age_screen # age at screens
    
    if(censor_type == "clinical"){ # if clinical, add tau_PC as endpoint of last interval
      age_screen <- c(age_screen, d_obs_censor_i$censor_time)
    }else if(censor_type == "censored"){ # if censored, add Inf as endpoint of last interval
      age_screen <- c(age_screen, Inf)
    }
    
    dt         <- diff(age_screen) # length of intervals
    K          <- length(dt)       # number of intervals
    
    # construct probability vector p
    prob <- if(censor_type == "screen"  ){
      
      theta$beta * (1-theta$beta)^((K-1):0) * dt[1:K]
      
    }  else if(censor_type == "clinical"){
      
      (1-theta$beta)^((K-1):0) * dt[1:K]
      
    } else if(censor_type == "censored"){
      
      prob_sojourn_H <- pweibull2_ab(
        age_screen[1:K    ] - 40, 
        age_screen[1:K + 1] - 40, 
        shape = theta$shape_H, rate = theta$rate_H, 
        log = FALSE
      )
      
      prob_screens <- (1-theta$beta)^((K-1):0)
      
      prob_sojourn_H * prob_screens
      
    }
    
    prob <- prob/sum(prob) # normalize probabilities
      
    } # end-same code as in rprop_tau_HP
  
  # find k_new
  k_new <- sum(age_screen < tau_HP) 
  
  if(censor_type == "censored" & k_new == K){
    log(prob[k_new]) # log(categorical(p)
  }else{
    log(prob[k_new]) - log(age_screen[k_new+1] - age_screen[k_new]) # log(categorical(p) * uniform(dt_k))
  }
  
}

MH_tau_PH <- function(d_obs_screen_i, d_obs_censor_i, indolent_i, tau_HP_i, theta){
  
  # data for person i
  #d_obs_screen_i <- d_obs_screen %>% filter(person_id == i)
  #d_obs_censor_i <- d_obs_censor %>% filter(person_id == i)
  #indolent_i     <- indolent[i] 
  
  # propose new tau_HP
  tau_HP_cur <- tau_HP_i
  tau_HP_new <- rprop_tau_HP(d_obs_screen_i, d_obs_censor_i, theta)
  if(tau_HP_new == tau_HP_cur)  return(tau_HP_cur)
  
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
MCMC <- function(m=100){
  
  {
  n_obs <- nrow(d_process) 
  
  RATE_H <- RATE_P <- BETA <- PSY <- numeric(m)
  TAU_HP <- matrix(nrow = m, ncol = n_obs)
  indolent <- d_process$indolent
  tau_HP_tmp <- tau_HP   <- if_else(d_process$tau_HP < d_obs_censor$censor_time, d_process$tau_HP, Inf)
  indolent               <- d_process$indolent #if_else(d_process$tau_HP < d_obs_censor$censor_time, d_process$indolent, NA_real_)
  theta_tmp <- theta
  
  n_screen_positive <- sum(d_obs_screen$screen_detected)
  
  # pre-process data
  d_id <- d_process %>% select(person_id)
  d_obs_screen_tbl <- d_obs_screen %>% nest(screens = screen_id:screen_detected)
  d_obs_censor_tbl <- d_obs_censor %>% nest(censor = censor_type:censor_time)
  screens   <- d_obs_screen_tbl$screens
  censoring <- d_obs_censor_tbl$censor 
  
  
  tic()
  for(n in 1 : m){ print(n)
    
    # update parameters
    rate_H    <- gibbs_rate_H(tau_HP, theta_tmp, prior, d_obs_censor)
    rate_P    <- gibbs_rate_P(tau_HP, indolent, theta_tmp, prior, d_obs_censor)
    psy       <- gibbs_psy   (indolent , prior)
    beta      <- gibbs_beta  (tau_HP, d_id, d_obs_screen, prior, n_screen_positive)
    theta_tmp <- list(rate_H=rate_H, rate_P=rate_P, psy=psy, beta=beta,
                      shape_H=theta_tmp$shape_H, shape_P=theta_tmp$shape_P)
    
    # update latent data
    #foreach(i=1:n_obs, .packages="tidyverse") %dopar% {
    #  source("helpers.R")
    tau_HP <- list(screens, censoring, indolent, tau_HP) %>%
      pmap_dbl(MH_tau_PH, theta_tmp)
    
    #for(i in 1 : n_obs){ # this for-loop should be done in parallel
    #  tau_HP_test2[i] <- MH_tau_PH(
    #    d_obs_screen_tbl$screens[[i]], d_obs_censor[i,], indolent[i], tau_HP[i], theta_tmp
    #    )
      #indolent[i] <- MH_indolent()
    #}
    
    #print(paste0("prop of Inf:", mean(is.infinite(tau_HP_tmp) & is.infinite(tau_HP)),
    #             " - prop of update:", mean(tau_HP_tmp != tau_HP)
    #             ))
    #tau_HP_tmp <- tau_HP
    
    RATE_H[n ] <- rate_H
    RATE_P[n ] <- rate_P
    BETA  [n ] <- beta
    PSY   [n ] <- psy
    TAU_HP[n,] <- tau_HP
    
  }
  
  runtime <- toc()
  } # end runtime
  
    # output
  THETA <- list(
    RATE_H=RATE_H, RATE_P=RATE_P, BETA=BETA, PSY=PSY
    )
  out <- list(
    THETA=THETA, TAU_HP=TAU_HP, runtime=runtime
    )
  return(out)
  
}

