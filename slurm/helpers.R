
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

rweibull_trunc <- function(a, b, shape, scale){
  u <- runif(1, pweibull(a, shape, scale), pweibull(b, shape, scale))
  qweibull(u, shape, scale)
}

dweibull_trunc <- function(x, a, b, shape, scale, log = T){
  if(log){
    dweibull(x, shape, scale, log = T) - log(pweibull_ab(a, b, shape, scale))
  }else{
    dweibull(x, shape, scale, log = F) /     pweibull_ab(a, b, shape, scale)
  }
}

#
# Biological process ####

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
# Gibbs theta ####

gibbs_rate_H <- function(tau_HP, theta, prior, censor_time){ # gamma prior and weibull likelihood conjugacy
  
  a_n <- prior$shape_H + sum(is.finite(tau_HP)) # number of observed weibull
  
  sojourn_H   <- pmin(tau_HP, censor_time) - 40
  sum_sojourn <- sum(sojourn_H^theta$shape_H) 
  b_n         <- prior$rate_H + sum_sojourn
  
  rgamma(1, shape = a_n, rate = b_n)
}

gibbs_rate_P <- function(tau_HP, indolent, theta, prior, censor_type, censor_time){ # gamma prior and weibull likelihood conjugacy
  
  a_n <- prior$shape_P + sum(censor_type == "clinical") # number of observed weibull
  
  sojourn_P   <- (censor_time - tau_HP)[!indolent & is.finite(tau_HP)]
  sum_sojourn <- sum(sojourn_P^theta$shape_P) 
  b_n         <- prior$rate_P + sum_sojourn
  
  rgamma(1, shape = a_n, rate = b_n)
}

gibbs_psi <- function(indolent, prior){ # beta prior and binomial likelihood conjugacy
  
  a_n <- prior$a_psi + sum(indolent  , na.rm = TRUE) 
  b_n <- prior$b_psi + sum(1-indolent, na.rm = TRUE) 
  
  rbeta(1, a_n, b_n)
}

gibbs_beta <- function(tau_HP, age_screen, prior, n_screen_positive){ # beta prior and binomial likelihood conjugacy
  
  n_screen <- list(tau_HP, age_screen) %>% 
    pmap_dbl(~sum(..1 < ..2)) %>%
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

dlog_psi <- function(indolent_i, theta){
  if(is.na(indolent_i)){
    0
  } else if(indolent_i) { # Bernoulli pmf
    log(theta$psi)
  } else if(!indolent_i){
    log(1-theta$psi)
  }
  
    # dbinom(
    #   indolent_i, 1,
    #   prob = theta$psi, log = TRUE
    # )
}

dlog_beta <- function(age_screen_i, theta, tau_HP_i, n_screen_positive_i){
  
  n_screen_i <- sum(age_screen_i > tau_HP_i) # number of screens during pre-clinical phase
  
  # product of Bernoulli RVs
  n_screen_positive_i*log(theta$beta) + (n_screen_i - n_screen_positive_i) * log(1-theta$beta)
  
}

dlog_likelihood_i <- function(tau_HP_i, indolent_i, censor_type_i, censor_time_i, age_screen_i, n_screen_positive_i, theta){
  
  dlog_H <- dlog_sojourn_H(tau_HP_i, theta, censor_time_i)
  dlog_P <- dlog_sojourn_P(tau_HP_i, theta, indolent_i, censor_type_i, censor_time_i)
  dlog_b <- dlog_beta     (age_screen_i, theta, tau_HP_i, n_screen_positive_i)
  dlog_p <- dlog_psi      (indolent_i, theta)
  
  dlog_likelihood <- dlog_H + dlog_P + dlog_b + dlog_p
  return(dlog_likelihood)
  
}

dlog_likelihood <- function(tau_HP, indolent, censor_type, censor_time, age_screen, n_screen_positive, theta, n_cpu){
  mcmapply(
    dlog_likelihood_i,
    tau_HP, indolent, censor_type, censor_time, age_screen, n_screen_positive,
    MoreArgs = list(theta = theta),
    mc.cores = n_cpu
  ) %>% sum()
}

dloglik_indolent_i <- function(tau_HP_i, indolent_i, censor_type_i, censor_time_i, theta){
  
  dlog_P <- dlog_sojourn_P(tau_HP_i, theta, indolent_i, censor_type_i, censor_time_i)
  dlog_p <- dlog_psi      (indolent_i, theta)
  
  dlog_likelihood <- dlog_P + dlog_p
  return(dlog_likelihood)
  
}

dloglik_indolent <- function(tau_HP, indolent, censor_type, censor_time, theta, n_cpu){
  mcmapply(
    dloglik_indolent_i,
    tau_HP, indolent, censor_type, censor_time,
    MoreArgs = list(theta = theta),
    mc.cores = n_cpu
  ) %>% sum()
}


#
# M-H tau_PH ####

compute_endpoints_i <- function(age_screen, censor_type_i, censor_time_i){
  
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

compute_endpoints <- function(age_screen, censor_type, censor_time, n_cpu){
  mcmapply(
    compute_endpoints_i,
    age_screen, censor_type, censor_time,
    mc.cores = n_cpu
  )
}

compute_prob_tau_i <- function(censor_type_i, censor_time_i, endpoints_i, theta){

  K <- length(endpoints_i) - 1 # number of intervals
  
  # compute prob of interval based on  (ignoring the mommograms)
  prob_interval <- if(censor_type_i == "screen"  ){
    
    pweibull_ab(
      endpoints_i[1:K    ] - 40,
      endpoints_i[1:K + 1] - 40,
      shape = theta$shape_H, scale = theta$scale_H
    )
    
  } else if(censor_type_i == "clinical"){
    
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
  
  prob_tau <- prob_interval * prob_screens
  prob_tau <- prob_tau/sum(prob_tau) # normalize probabilities
  
  return(prob_tau)
  
}

compute_prob_tau <- function(censor_type, censor_time, endpoints, theta, n_cpu){
  mcmapply(
    compute_prob_tau_i,
    censor_type, censor_time, endpoints,
    MoreArgs = list(theta = theta),
    USE.NAMES = FALSE,
    mc.cores = n_cpu
  )
}

rprop_tau_HP_i <- function(censor_type_i, censor_time_i, endpoints_i, prob_tau_i, theta){
  
  K <- length(prob_tau_i) # number of intervals
  
  # sample interval
  k_new <- sample.int(K, 1, prob = prob_tau_i) # sample interval
  
  # sample tau_HP in chosen interval
  if(censor_type_i == "screen"  ){
    
    sojourn_H_new <- rweibull_trunc(
      endpoints_i[k_new] - 40, endpoints_i[k_new+1] - 40,
      theta$shape_H, theta$scale_H
      )
    tau_HP_new <- 40 + sojourn_H_new
    
  } else if(censor_type_i == "clinical"){
    
    sojourn_P_new <- rweibull_trunc(
      censor_time_i - endpoints_i[k_new+1], censor_time_i - endpoints_i[k_new],
      theta$shape_P, theta$scale_P
      )
    tau_HP_new <- censor_time_i - sojourn_P_new
    
  } else if(censor_type_i == "censored"){
    
    if(k_new == K){ # tau_HP_new > time_censor
      
      tau_HP_new <- Inf
      
    }else{ # tau_HP_new < time_censor
      
      sojourn_H_new <- rweibull_trunc(
        endpoints_i[k_new] - 40, endpoints_i[k_new+1] - 40,
        theta$shape_H, theta$scale_H
      )
      tau_HP_new <- 40 + sojourn_H_new
      
    }
    
  }
  
  return(tau_HP_new)

}

rprop_tau_HP <- function(censor_type, censor_time, endpoints, prob_tau, theta, n_cpu){
  mcmapply(
    rprop_tau_HP_i,
    censor_type, censor_time, endpoints, prob_tau,
    MoreArgs = list(theta = theta),
    USE.NAMES = FALSE,
    mc.cores = n_cpu
  )
}

dlog_prop_tau_HP_i <- function(tau_HP_i, censor_type_i, censor_time_i, endpoints_i, prob_tau_i, theta){
  
  K      <- length(prob_tau_i) # number of intervals
  
  # contribution of k_new
  k_new  <- sum(endpoints_i < tau_HP_i)
  dlog_k <- log(prob_tau_i[k_new])
  
  # contribution of tau_HP_i
  dlog_tau <- if(censor_type_i == "screen"  ){
    
    sojourn_H <- tau_HP_i - 40
    dweibull_trunc(
      sojourn_H, endpoints_i[k_new] - 40, endpoints_i[k_new+1] - 40,
      theta$shape_H, theta$scale_H, log = T
      )
    
  } else if(censor_type_i == "clinical"){
    
    sojourn_P <- censor_time_i - tau_HP_i
    dweibull_trunc(
      sojourn_P, censor_time_i - endpoints_i[k_new+1], censor_time_i - endpoints_i[k_new],
      theta$shape_P, theta$scale_P, log = T
    )
    
  } else if(censor_type_i == "censored"){
    
    if(k_new == K){ # tau_HP_i > time_censor
      
      0
      
    }else{ # tau_HP_i < time_censor
      
      sojourn_H <- tau_HP_i - 40
      dweibull_trunc(
        sojourn_H, endpoints_i[k_new] - 40, endpoints_i[k_new+1] - 40,
        theta$shape_H, theta$scale_H, log = T
      )
      
    }
    
  }
  
  return(dlog_k + dlog_tau)
  
}

#
# indolent ####

compute_prob_indolent_i <- function(tau_HP_i, censor_type_i, censor_time_i, theta){
  
  L_0 <- dloglik_indolent_i(tau_HP_i, indolent_i = 0, censor_type_i, censor_time_i, theta) %>% exp()
  L_1 <- dloglik_indolent_i(tau_HP_i, indolent_i = 1, censor_type_i, censor_time_i, theta) %>% exp()
  
  p_i <- L_1 / (L_0 + L_1)
  
}

compute_prob_indolent <- function(tau_HP, censor_type, censor_time, theta, n_cpu){
  mcmapply(
    compute_prob_indolent_i,
    tau_HP, censor_type, censor_time,
    MoreArgs = list(theta = theta),
    USE.NAMES = FALSE,
    mc.cores = n_cpu
  )
}

rprop_indolent_i <- function(tau_HP_i, censor_type_i, prob_indolent_i){
  
  indolent_i <- if(censor_type_i == "clinical"){
    0
  }else if(is.infinite(tau_HP_i)){
    NA
  }else{
    #prob_indolent_i <- compute_prob_indolent(tau_HP_i, censor_type_i, censor_time_i, theta)
    rbinom(1,1,prob_indolent_i)
  }
  
  return(indolent_i)
  
}

rprop_indolent <- function(tau_HP, censor_type, prob_indolent, n_cpu){
  mcmapply(
    rprop_indolent_i,
    tau_HP, censor_type, prob_indolent,
    #MoreArgs = list(theta = theta),
    USE.NAMES = FALSE,
    mc.cores = n_cpu
  )
}


dlog_prop_indolent_i <- function(indolent_i, tau_HP_i, censor_type_i, prob_indolent_i){
  
  dlog <- if(censor_type_i == "clinical"){
    0
  }else if(is.infinite(tau_HP_i)){
    0
  }else if(indolent_i){ # Bernoulli pmf
    log(prob_indolent_i)
  }else if(!indolent_i){
    log(1-prob_indolent_i)
  }
  
    #prob_indolent_i <- compute_prob_indolent(tau_HP_i, censor_type_i, censor_time_i, theta)
    #dbinom(indolent_i, 1, prob_indolent_i, log = TRUE)
  
  return(dlog)
  
}

dlog_prop_indolent <- function(indolent, tau_HP, censor_type, prob_indolent, n_cpu){
  mcmapply(
    dlog_prop_indolent_i,
    indolent, tau_HP, censor_type, prob_indolent,
    #MoreArgs = list(theta = theta),
    mc.cores = n_cpu
  ) %>% sum()
}


#
# psi ####
rprop_psi <- function(theta_cur, epsilon){
  
  psi_prop <- runif(1, theta_cur$psi-epsilon, theta_cur$psi+epsilon)
  
  # reflection on lower bound 0 and upper bound 1 
  if(psi_prop<0){
    psi_prop <- 0 + (0-psi_prop)
  }else if(psi_prop>1){
    psi_prop <- 1 - (psi_prop-1)
  } 
  
  return(psi_prop)
  
}

# latent ####

rprop_latent_i <- function(censor_type_i, censor_time_i, endpoints_i, prob_tau_i, prob_indolent_i, theta){

  # propose tau
  tau_HP_new <- rprop_tau_HP_i(censor_type_i, censor_time_i, endpoints_i, prob_tau_i, theta)
  
  # propose indolent
  indolent_new <- rprop_indolent_i(tau_HP_new, censor_type_i, prob_indolent_i)
  
  #latent_new <- c("tau_HP" = tau_HP_new, "indolent" = indolent_new)
  latent_new <- c(tau_HP_new, indolent_new)
  
  return(latent_new)
  
}

rprop_latent <- function(
    censor_type, censor_time, endpoints, prob_tau, prob_indolent, theta, n_cpu
    ){
  mcmapply(
    rprop_latent_i,
    censor_type, censor_time, endpoints, prob_tau, prob_indolent,
    MoreArgs = list(theta = theta),
    USE.NAMES = FALSE,
    mc.cores = n_cpu
  )
}
  
dlog_prop_latent_i <- function(
    tau_HP_i, indolent_i, censor_type_i, censor_time_i, endpoints_i,
    prob_tau_i, prob_indolent_i, theta
    ){
  
  # dlog tau
  dlog_tau <- dlog_prop_tau_HP_i(tau_HP_i, censor_type_i, censor_time_i, endpoints_i, prob_tau_i, theta)
  
  # dlog indolent
  dlog_indolent <- dlog_prop_indolent_i(indolent_i, tau_HP_i, censor_type_i, prob_indolent_i)
  
  dlog_prop <- dlog_tau + dlog_indolent
  return(dlog_prop)
  
}

dlog_prop_latent <- function(tau_HP, indolent, censor_type, censor_time, endpoints, prob_tau, prob_indolent, theta, n_cpu){
  mcmapply(
    dlog_prop_latent_i,
    tau_HP, indolent, censor_type, censor_time, endpoints, prob_tau, prob_indolent,
    MoreArgs = list(theta = theta),
    USE.NAMES = FALSE,
    mc.cores = n_cpu
  ) %>% sum()
}


#
# M-H (tau, indolent) ####

MH_latent_i <- function(tau_HP_cur_i, indolent_cur_i, censor_type_i, censor_time_i, age_screen_i, endpoints_i, n_screen_positive_i, theta){
  
  # propose new latent
  prob_tau_i          <- compute_prob_tau_i(censor_type_i, censor_time_i, endpoints_i, theta)
  tau_HP_new_i        <- rprop_tau_HP_i(censor_type_i, censor_time_i, endpoints_i, prob_tau_i, theta)
  prob_indolent_new_i <- compute_prob_indolent_i(tau_HP_new_i, censor_type_i, censor_time_i, theta)
  indolent_new_i      <- rprop_indolent_i(tau_HP_new_i, censor_type_i, prob_indolent_new_i)
  
  # M-H acceptance ratio
  prob_indolent_cur_i <- compute_prob_indolent_i(tau_HP_cur_i, censor_type_i, censor_time_i, theta)
  dlog_prop_cur_i     <- dlog_prop_latent_i(tau_HP_cur_i, indolent_cur_i, censor_type_i, censor_time_i, endpoints_i, prob_tau_i, prob_indolent_cur_i, theta)
  dlog_prop_new_i     <- dlog_prop_latent_i(tau_HP_new_i, indolent_new_i, censor_type_i, censor_time_i, endpoints_i, prob_tau_i, prob_indolent_new_i, theta)
  dlog_lik_cur_i      <- dlog_likelihood_i(tau_HP_cur_i, indolent_cur_i, censor_type_i, censor_time_i, age_screen_i, n_screen_positive_i, theta)
  dlog_lik_new_i      <- dlog_likelihood_i(tau_HP_new_i, indolent_new_i, censor_type_i, censor_time_i, age_screen_i, n_screen_positive_i, theta)
  
  MH_logratio_i <- dlog_lik_new_i - dlog_lik_cur_i + dlog_prop_cur_i - dlog_prop_new_i
  
  out <- if(runif(1) < exp(MH_logratio_i)){ # accept new latent data
    c(tau_HP_new_i, indolent_new_i, TRUE)
  }else{ # keep current latent data
    c(tau_HP_cur_i, indolent_cur_i, FALSE)
  }
  
  return(out)
  
}

MH_latent <- function(tau_HP_cur, indolent_cur, censor_type, censor_time, age_screen, endpoints, n_screen_positive, theta, n_cpu){
  mcmapply(
    MH_latent_i,
    tau_HP_cur, indolent_cur, censor_type, censor_time, age_screen, endpoints, n_screen_positive,
    MoreArgs = list(theta = theta),
    USE.NAMES = FALSE,
    mc.cores = n_cpu
  )
}

#
# M-H (psi, indolent) ####

rprop_dlog_indolent_i <- function(tau_HP_i, censor_type_i, censor_time_i, theta){
  
  prob_indolent_new_i <- compute_prob_indolent_i(tau_HP_i, censor_type_i, censor_time_i, theta)
  indolent_new_i      <- rprop_indolent_i(tau_HP_i, censor_type_i, prob_indolent_new_i)
  dlog_prop_new_i     <- dlog_prop_indolent_i(indolent_new_i, tau_HP_i, censor_type_i, prob_indolent_new_i)
  dlog_lik_new_i      <- dloglik_indolent_i(tau_HP_i, indolent_new_i, censor_type_i, censor_time_i, theta)
  
  out <- c(indolent_new_i, dlog_prop_new_i, dlog_lik_new_i)
  return(out)
  
}

rprop_dlog_indolent <- function(tau_HP, censor_type, censor_time, theta, n_cpu){
  
  out_mc <- mcmapply(
    rprop_dlog_indolent_i,
    tau_HP, censor_type, censor_time,
    MoreArgs = list(theta = theta),
    USE.NAMES = FALSE,
    mc.cores = n_cpu
  )
  
  indolent_new           <- out_mc[1,]
  dlog_prop_indolent_new <- out_mc[2,] %>% sum()
  dlog_lik_new           <- out_mc[3,] %>% sum()
  
  out <- list(
    indolent_new=indolent_new, 
    dlog_prop_indolent_new=dlog_prop_indolent_new, 
    dlog_lik_new=dlog_lik_new)
  return(out)
  
}



MH_psi_indolent <- function(
    theta_cur, indolent_cur, tau_HP, censor_type, censor_time, epsilon, prior, n_cpu
    ){
  
  # propose psi
  psi_new <- rprop_psi(theta_cur, epsilon) # symmetric proposal
  theta_new <- theta_cur
  theta_new[["psi"]] <- psi_new
  
  # prior density
  dlog_prior_cur <- dbeta(theta_cur$psi, prior$a_psi, prior$b_psi, log=T)
  dlog_prior_new <- dbeta(theta_new$psi, prior$a_psi, prior$b_psi, log=T)
  
  # propose indolent, compute log likelihood and log proposal density
  out <- rprop_dlog_indolent(tau_HP, censor_type, censor_time, theta_new, n_cpu)
  #indolent_new           <- rprop_indolent(tau_HP, censor_type, prob_indolent_new, n_cpu)
  
  # proposal density
  prob_indolent_cur      <- compute_prob_indolent(tau_HP, censor_type, censor_time, theta_cur, n_cpu)
  dlog_prop_indolent_cur <- dlog_prop_indolent(indolent_cur, tau_HP, censor_type, prob_indolent_cur, n_cpu)
  dlog_prop_indolent_new <- out[["dlog_prop_indolent_new"]]
  #prob_indolent_new      <- compute_prob_indolent(tau_HP, censor_type, censor_time, theta_new, n_cpu)
  #dlog_prop_indolent_new <- dlog_prop_indolent(indolent_new, tau_HP, censor_type, prob_indolent_new, n_cpu)
  
  # log likelihood
  dlog_lik_cur <- dloglik_indolent(tau_HP, indolent_cur, censor_type, censor_time, theta_cur, n_cpu)
  dlog_lik_new <- out[["dlog_lik_new"]]
  #dlog_lik_new <- dloglik_indolent(tau_HP, indolent_new, censor_type, censor_time, theta_new, n_cpu)

  # M-H acceptance ratio
  MH_logratio <- dlog_lik_new   - dlog_lik_cur + 
                 dlog_prior_new - dlog_prior_cur + 
                 dlog_prop_indolent_cur  - dlog_prop_indolent_new 
  
  out <- if(runif(1) < exp(MH_logratio)){ # accept new values
    list(indolent = out[["indolent_new"]], theta = theta_new, accept = TRUE)
  } else { # keep current values
    list(indolent = indolent_cur, theta = theta_cur, accept = FALSE)
  }

  return(out)
  
}

#
# MCMC ####

MCMC <- function(
    d_process, d_obs_screen, d_obs_censor,
    theta_0, prior, epsilon,
    M = 100, thin = 1,
    n_cpu = 1, verbose = FALSE
    ){
  
  {
    
  # setup
  n_obs <- nrow(d_process)
  
  M_thin <- M / thin
  RATE_H <- RATE_P <- BETA <- PSI <- ACCEPT_PSI <- numeric(M_thin)
  TAU_HP <- INDOLENT <- ACCEPT_LATENT <- matrix(nrow = M_thin, ncol = n_obs)
  
  # pre-process data
  theta <- theta_0  
  censor_type             <- d_obs_censor$censor_type
  n_screen_positive       <- censor_type=="screen"
  n_screen_positive_total <- sum(d_obs_screen$screen_detected)
  censor_time             <- d_obs_censor$censor_time
  d_obs_screen_tbl        <- d_obs_screen %>% nest(screens = screen_id:screen_detected)
  screens                 <- d_obs_screen_tbl$screens
  age_screen              <- screens %>% map(~ .[["age_screen"]])
  endpoints               <- compute_endpoints(age_screen, censor_type, censor_time, n_cpu)

  # initialisation of latent variables
  prob_tau_0      <- compute_prob_tau(censor_type, censor_time, endpoints, theta, n_cpu)
  prob_indolent_0 <- compute_prob_indolent(tau_HP, censor_type, censor_time, theta, n_cpu)
  latent_data <- rprop_latent(censor_type, censor_time, endpoints, prob_tau_0, prob_indolent_0, theta, n_cpu)
  tau_HP      <- latent_data[1,]
  indolent    <- latent_data[2,]
  
  # MCMC
  tic()
  for(m in 1 : M){ 
    if(verbose)  print(paste0("n=", n_obs, "-M=", M, "-n_cpu=", n_cpu, "-epsilon=", epsilon, " -- ", m, "/", M, " = ", round(100*m/M, 2), "%"))
    
    # update (rate_H, rate_P, beta)
    rate_H <- gibbs_rate_H(tau_HP, theta, prior, censor_time)
    rate_P <- gibbs_rate_P(tau_HP, indolent, theta, prior, censor_type, censor_time)
    beta   <- gibbs_beta  (tau_HP, screens, prior, n_screen_positive_total)
    theta  <- list(
      rate_H  = rate_H         , rate_P  = rate_P         , psi = theta$psi, beta = beta,
      shape_H = theta_0$shape_H, shape_P = theta_0$shape_P
      ) %>% update_scales()
    
    # update (psi, indolent)
    out_psi_indolent <- MH_psi_indolent(theta, indolent, tau_HP, censor_type, censor_time, epsilon, prior, n_cpu)
    indolent         <- out_psi_indolent[["indolent"]]
    theta            <- out_psi_indolent[["theta"]]
    accept_psi       <- out_psi_indolent[["accept"]]
    
    # update (tau_HP, indolent)
    out_latent    <- MH_latent(tau_HP, indolent, censor_type, censor_time, age_screen, endpoints, n_screen_positive, theta, n_cpu)
    tau_HP        <- out_latent[1,]
    indolent      <- out_latent[2,]
    accept_latent <- out_latent[3,]
    
    # save output
    if(m %% thin == 0){
      
      m_thin <- m %/% thin
      
      RATE_H  [m_thin ] <- rate_H
      RATE_P  [m_thin ] <- rate_P
      BETA    [m_thin ] <- beta
      PSI     [m_thin ] <- theta$psi
      TAU_HP  [m_thin,] <- tau_HP
      INDOLENT[m_thin,] <- indolent
      ACCEPT_LATENT[m_thin,] <- accept_latent
      ACCEPT_PSI   [m_thin ] <- accept_psi
      
    }
    
  }
  
  tictoc  <- toc()
  runtime <- tictoc$toc - tictoc$tic
  } # end runtime
  
  # output
  THETA <- list(
    RATE_H=RATE_H, RATE_P=RATE_P, BETA=BETA, PSI=PSI
    )
  
  out <- list(
    THETA=THETA, TAU_HP=TAU_HP, INDOLENT=INDOLENT,
    ACCEPT_LATENT=ACCEPT_LATENT, 
    ACCEPT_PSI=ACCEPT_PSI,
    epsilon=epsilon,
    runtime=runtime
    )
  return(out)
  
}
