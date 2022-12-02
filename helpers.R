
#
# Weibull ####

rate2scale <- function(rate, shape){
  rate^(-1/shape)
}

update_scales <- function(theta){
  theta$scale_H <- rate2scale(theta$rate_H, theta$shape_H)
  theta$scale_P <- rate2scale(theta$rate_P, theta$shape_P)
  return(theta)
}

pweibull_ab <- function(a, b, shape, scale){
  pweibull(b, shape, scale, lower.tail = T) - 
  pweibull(a, shape, scale, lower.tail = T)
}

rweibul_trunc <- function(a, b, shape, scale){
  u <- runif(1, pweibull(a, shape, scale), pweibull(b, shape, scale))
  qweibull(u, shape, scale)
}

dweibul_trunc <- function(x, a, b, shape, scale, log = T){
  if(log){
    dweibull(x, shape, scale, log = T) - log(pweibull_ab(a, b, shape, scale))
  }else{
    dweibull(x, shape, scale, log = F) /     pweibull_ab(a, b, shape, scale)
  }
}

#
# Process ####

draw_sojourn_H <- function(theta){
  rweibull(1, shape = theta$shape_H, scale = theta$scale_H)
}

draw_sojourn_P <- function(theta){
  rweibull(1, shape = theta$shape_P, scale = theta$scale_P)
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

gibbs_rate_H <- function(tau_HP, theta, prior, censor_time){ # gamma prior and weibull likelihood conjugacy
  
  a_n <- prior$shape_H + sum(is.finite(tau_HP)) # number of observed weibull
  
  sojourn_H   <- pmin(tau_HP, censor_time) - 40
  sum_sojourn <- sum(sojourn_H^theta$shape_H) 
  b_n         <- prior$rate_H + sum_sojourn
  
  rgamma(1, shape = a_n, rate = b_n)
}

gibbs_rate_P <- function(tau_HP, indolent, theta, prior, censor_type, censor_time){ # gamma prior and weibull likelihood conjugacy
  
  a_n <- prior$shape_P + sum(censor_type == "clinical") # numbe rof observed weibull
  
  sojourn_P   <- (censor_time - tau_HP)[!indolent & is.finite(tau_HP)]
  sum_sojourn <- sum(sojourn_P^theta$shape_P) 
  b_n         <- prior$rate_P + sum_sojourn
  
  rgamma(1, shape = a_n, rate = b_n)
}

gibbs_psy <- function(indolent, prior){ # beta prior and binomial likelihood conjugacy
  
  a_n <- prior$a_psy + sum(indolent  , na.rm = TRUE) 
  b_n <- prior$b_psy + sum(1-indolent, na.rm = TRUE) 
  
  rbeta(1, a_n, b_n)
}

gibbs_beta <- function(tau_HP, screens, prior, n_screen_positive){ # beta prior and binomial likelihood conjugacy
  
  n_screen <- list(tau_HP, screens) %>% 
    pmap_dbl(~sum(..1 < ..2$age_screen)) %>%
    sum()
  
  a_n <- prior$a_beta + n_screen_positive
  b_n <- prior$b_beta + (n_screen-n_screen_positive)
  
  rbeta(1, a_n, b_n)
}

#
# Likelihood ####

dlog_sojourn_H <- function(tau_HP_i, theta, censor_time_i){
  
  if(is.infinite(tau_HP_i)){
    pweibull(
      censor_time_i - 40,
      shape = theta$shape_H, scale = theta$scale_H,
      log.p = TRUE, lower.tail = FALSE
      )
  }else if(is.finite(tau_HP_i)){
    dweibull(
      tau_HP_i - 40, 
      shape = theta$shape_H, scale = theta$scale_H,
      log = TRUE
    )
  }
  
}

dlog_sojourn_P <- function(tau_HP_i, theta, indolent_i, censor_type_i, censor_time_i){
  
  if(is.infinite(tau_HP_i)){
    0
  }else if(indolent_i){
    0
  }else if(censor_type_i == "clinical"){
    dweibull(
      censor_time_i - tau_HP_i, 
      shape = theta$shape_P, scale = theta$scale_P, 
      log = TRUE
    )
  } else if(censor_type_i %in% c("screen", "censored")){
    pweibull(
      censor_time_i - tau_HP_i, 
      shape = theta$shape_P, scale = theta$scale_P, 
      log.p = TRUE, lower.tail = FALSE
    )
  }
  
}

dlog_psy <- function(indolent_i, theta){
  if(is.na(indolent_i)){
    0
  } else {
    dbinom(
      indolent_i, 1, 
      prob = theta$psy, log = TRUE
    )
  }
}

dlog_beta <- function(d_obs_screen_i, theta, tau_HP_i, censor_type_i){
  
  n_screen          <- sum(d_obs_screen_i$age_screen>tau_HP_i) # number of screens during pre-clinical  
  n_screen_positive <- censor_type_i=="screen" # number of positive screens during pre-clinical  
  
  # product of Bernoulli RVs
  n_screen_positive*log(theta$beta) + (n_screen - n_screen_positive)*log(1-theta$beta)
  
}

dlog_likelihood_i <- function(tau_HP_i, indolent_i, d_obs_screen_i, censor_type_i, censor_time_i, theta){
  
  dlog_H <- dlog_sojourn_H(tau_HP_i, theta, censor_time_i)
  dlog_P <- dlog_sojourn_P(tau_HP_i, theta, indolent_i, censor_type_i, censor_time_i)
  dlog_p <- dlog_psy      (indolent_i, theta)
  dlog_b <- dlog_beta     (d_obs_screen_i, theta, tau_HP_i, censor_type_i)
  
  dlog_likelihood <- dlog_H + dlog_P + dlog_p + dlog_b
  return(dlog_likelihood)
  
}


#
# MH ####

compute_endpoints <- function(screens_i, censor_type_i, censor_time_i){
  
  age_screen <- screens_i$age_screen # age at screens
  
  # construct endpoints of intervals (t_i, t_{i+1}]
  endpoints <- if(censor_type_i == "screen"){ # if screen, age_screen are endpoints
    
    age_screen 
    
  }else if(censor_type_i == "clinical"){ # if clinical, add tau_PC as endpoint of last interval
    
    c(age_screen, censor_time_i     )
    
  }else if(censor_type_i == "censored"){ # if censored, add tau_PC and Inf as endpoints
    
    c(age_screen, censor_time_i, Inf)
    
  }
  
  return(endpoints)
  
}

compute_prob <- function(censor_type_i, censor_time_i, endpoints_i, theta){
  
  
  K <- length(endpoints_i) - 1 # number of intervals
  
  # compute prob of interval based on  (ignoring the mommograms)
  prob_interval <- if(censor_type_i == "screen"  ){
    
    #diff(endpoints_i) # length of intervals
    pweibull_ab(
      endpoints_i[1:K    ] - 40,
      endpoints_i[1:K + 1] - 40,
      shape = theta$shape_H, scale = theta$scale_H
    )
    
  } else if(censor_type_i == "clinical"){
    
    #diff(endpoints_i) # length of intervals
    # pweibull_ab(
    #   endpoints_i[1:K    ] - 40,
    #   endpoints_i[1:K + 1] - 40,
    #   shape = theta$shape_H, scale = theta$scale_H
    # )
    pweibull_ab(
      censor_time_i - endpoints_i[1:K + 1],
      censor_time_i - endpoints_i[1:K    ],
      shape = theta$shape_P, scale = theta$scale_P
    )
    
  } else if(censor_type_i == "censored"){
    
    pweibull_ab(
      endpoints_i[1:K    ] - 40, 
      endpoints_i[1:K + 1] - 40, 
      shape = theta$shape_H, scale = theta$scale_H
    )
    
  }
  
  prob_screens <- if(censor_type_i == "screen"  ){
    
    (1-theta$beta)^((K-1):0) * theta$beta
    
  } else if(censor_type_i == "clinical"){
    
    (1-theta$beta)^((K-1):0)
    
  } else if(censor_type_i == "censored"){
    
    (1-theta$beta)^((K-1):0)
    
  }
  
  prob <- prob_interval * prob_screens
  prob <- prob/sum(prob) # normalize probabilities
  
  return(prob)
  
}

rprop_tau_HP <- function(censor_type_i, censor_time_i, endpoints_i, prob, theta){
  
  K <- length(prob) # number of intervals
  
  # sample interval
  k_new <- sample.int(K, 1, prob = prob) # sample interval
  
  # sample tau_HP in chosen interval
  if(censor_type_i == "screen"  ){
    
    sojourn_H_new <- rweibul_trunc(
      endpoints_i[k_new] - 40, endpoints_i[k_new+1] - 40,
      theta$shape_H, theta$scale_H
      )
    tau_HP_new <- 40 + sojourn_H_new
    
  } else if(censor_type_i == "clinical"){
    
    sojourn_P_new <- rweibul_trunc(
      censor_time_i - endpoints_i[k_new+1], censor_time_i - endpoints_i[k_new],
      theta$shape_P, theta$scale_P
      )
    tau_HP_new <- censor_time_i - sojourn_P_new
    
  } else if(censor_type_i == "censored"){
    
    if(k_new == K){
      
      tau_HP_new <- Inf # tau_HP > time_censor
      
    }else{ # tau_HP < time_censor
      
      sojourn_H_new <- rweibul_trunc(
        endpoints_i[k_new] - 40, endpoints_i[k_new+1] - 40,
        theta$shape_H, theta$scale_H
      )
      tau_HP_new <- 40 + sojourn_H_new
      
    }
    
  }
  
  return(tau_HP_new)

}

dlog_prop_tau_HP <- function(tau_HP, censor_type_i, censor_time_i, endpoints_i, prob, theta){
  
  K      <- length(prob) # number of intervals
  
  # contribution of k_new
  k_new  <- sum(endpoints_i < tau_HP)
  dlog_k <- log(prob[k_new])
  
  # contribution of tau_HP
  dlog_tau <- if(censor_type_i == "screen"  ){
    
    sojourn_H <- tau_HP - 40
    dweibul_trunc(
      sojourn_H, endpoints_i[k_new] - 40, endpoints_i[k_new+1] - 40,
      theta$shape_H, theta$scale_H, log = T
      )
    
  } else if(censor_type_i == "clinical"){
    
    sojourn_P <- censor_time_i - tau_HP
    dweibul_trunc(
      sojourn_P, censor_time_i - endpoints_i[k_new+1], censor_time_i - endpoints_i[k_new],
      theta$shape_P, theta$scale_P, log = T
    )
    
  } else if(censor_type_i == "censored"){
    
    if(k_new == K){
      
      0
      
    }else{ # tau_HP < time_censor
      
      sojourn_H <- tau_HP - 40
      dweibul_trunc(
        sojourn_H, endpoints_i[k_new] - 40, endpoints_i[k_new+1] - 40,
        theta$shape_H, theta$scale_H, log = T
      )
      
    }
    
  }
  
  return(dlog_k + dlog_tau)
  
}

MH_tau_PH <- function(d_obs_screen_i, censor_type_i, censor_time_i, endpoints_i, indolent_i, tau_HP_i, theta){
  
  # data for person i
  #d_obs_screen_i <- d_obs_screen %>% filter(person_id == i)
  #d_obs_censor_i <- d_obs_censor %>% filter(person_id == i)
  #indolent_i     <- indolent[i] 
  #tau_HP_i <- tau_HP[i]
  
  # propose new tau_HP
  tau_HP_cur <- tau_HP_i
  prob       <- compute_prob(censor_type_i, censor_time_i, endpoints_i, theta)
  tau_HP_new <- rprop_tau_HP(censor_type_i, censor_time_i, endpoints_i, prob, theta)
  if(tau_HP_new == tau_HP_cur)  return(tau_HP_cur)
  
  # M-H acceptance ratio
  dlog_prop_cur  <- dlog_prop_tau_HP (tau_HP_cur, censor_type_i, censor_time_i, endpoints_i, prob, theta)
  dlog_prop_new  <- dlog_prop_tau_HP (tau_HP_new, censor_type_i, censor_time_i, endpoints_i, prob, theta)
  dlog_lik_cur   <- dlog_likelihood_i(tau_HP_cur, indolent_i, d_obs_screen_i, censor_type_i, censor_time_i, theta)
  dlog_lik_new   <- dlog_likelihood_i(tau_HP_new, indolent_i, d_obs_screen_i, censor_type_i, censor_time_i, theta)

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
  
  
  # pre-process data
  n_screen_positive <- sum(d_obs_screen$screen_detected)
  d_id <- d_process %>% select(person_id)
  d_obs_screen_tbl <- d_obs_screen %>% nest(screens = screen_id:screen_detected)
  screens          <- d_obs_screen_tbl$screens
  censor_type      <- d_obs_censor$censor_type
  censor_time      <- d_obs_censor$censor_time
  endpoints        <- list(screens, censor_type, censor_time) %>% 
    pmap(compute_endpoints)
  
  # MCMC
  tic()
  for(n in 1 : m){ print(n)
    
    # update parameters
    rate_H    <- gibbs_rate_H(tau_HP, theta_tmp, prior, censor_time)
    rate_P    <- gibbs_rate_P(tau_HP, indolent, theta_tmp, prior, censor_type, censor_time)
    psy       <- gibbs_psy   (indolent , prior)
    beta      <- gibbs_beta  (tau_HP, screens, prior, n_screen_positive)
    theta_tmp <- list(rate_H=rate_H, rate_P=rate_P, psy=psy, beta=beta,
                      shape_H=theta_tmp$shape_H, shape_P=theta_tmp$shape_P)
    theta_tmp <- update_scales(theta_tmp)
    
    # update latent data
    #foreach(i=1:n_obs, .packages="tidyverse") %dopar% {
    tau_HP <- list(screens, censor_type, censor_time, endpoints, indolent, tau_HP) %>%
      pmap_dbl(MH_tau_PH, theta_tmp)
    
    # for(i in 1 : n_obs){ # this for-loop should be done in parallel
    #   set.seed(i)
    #   tau_HP_tmp[i] <- MH_tau_PH(
    #     screens[[i]], censor_type[i], censor_time[i], endpoints[[i]], indolent[i], tau_HP[i], theta_tmp
    #   )
    # }
    
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

