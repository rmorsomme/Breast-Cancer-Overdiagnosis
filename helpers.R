
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
    t < tau_PC ~ "P",
    TRUE       ~ "C"
  )
}

screen_result <- function(compartment, theta) {
  if(compartment == "H")  return(FALSE)  
  if(compartment == "P")  return(rbernoulli(1, p = theta$beta))
  if(compartment == "C")  print("error - screen when clinical"); return(NA)
}

#
# Gibbs ####

gibbs_rate_H <- function(d_process, theta, prior){
  
  a_0 <- prior$shape_H
  n   <- nrow(d_process)
  a_n <- a_0 + n
  
  b_0 <- prior$rate_H
  sojourn_H <- d_process$sojourn_H
  sum_sojourn <- sum(sojourn_H^theta$shape_H) 
  b_n <- b_0 + sum_sojourn
  
  rgamma(1, shape = a_n, rate = b_n)
}

gibbs_rate_P <- function(d_process, theta, prior){
  
  sojourn_P <- d_process %>% filter(indolent == 0) %>% pull(sojourn_P)
  
  a_0 <- prior$shape_P
  n   <- length(sojourn_P)
  a_n <- a_0 + n
  
  b_0 <- prior$rate_P
  sum_sojourn <- sum(sojourn_P^theta$shape_P) 
  b_n <- b_0 + sum_sojourn
  
  rgamma(1, shape = a_n, rate = b_n)
}

gibbs_psy <- function(d_process, prior){
  
  indolent <- d_process$indolent
  
  a_0 <- prior$a_psy
  x   <- sum(indolent) 
  a_n <- a_0 + x
  
  b_0 <- prior$b_psy
  n   <- length(indolent)
  b_n <- b_0 + (n-x)
  
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

dlog_sojourn_P <- function(sojourn_P, theta){
  dweibull2(
    sojourn_P, 
    shape = theta$shape_P, rate = theta$rate_P, 
    log = TRUE
  )
}

dlog_psy <- function(indolent, theta){
  dbinom(
    indolent, 1, 
    prob = theta$psy, log = TRUE
  )
}

dlog_beta <- function(d_obs_screen_i, theta, tau_HP, tau_PC){
  
  # Extract outcome of screens during preclinical phase
  screen_preclinical <- d_obs_screen_i %>%
    mutate(compartment = compute_compartment(age_screen, tau_HP, tau_PC)) %>%
    filter(compartment == "P") %>%
    pull(screen_detected)
  
  n_screen          <- length(screen_preclinical)
  n_screen_positive <- sum(screen_preclinical)
  
  # product of Bernoulli RVs
  n_screen_positive*log(theta$beta) + (n_screen - n_screen_positive)*log(1-theta$beta)
  
}

dlog_likelihood_i <- function(d_process_i, d_obs_screen_i, theta){
  
  dlog_H <- dlog_sojourn_H(d_process_i$sojourn_H, theta)
  dlog_P <- dlog_sojourn_P(d_process_i$sojourn_P, theta)
  dlog_p <- dlog_psy(d_process_i$indolent, theta)
  dlog_b <- dlog_beta(d_obs_screen_i, theta, d_process_i$tau_HP, d_process_i$tau_PC)
  
  dlog_likelihood <- dlog_H + dlog_P + dlog_p + dlog_b
  return(dlog_likelihood)
  
}


#
# MH ####
rprop_tau_PH <- function(d_obs_screen_i, d_obs_clinical_i){
  
  age_screen_positive <- d_obs_screen_i %>% filter(screen_detected == 1) %>% pull(age_screen)
  upper_bound         <- min(age_screen_positive, d_obs_clinical_i$tau_PC_obs, na.rm = TRUE)
  tau_PH              <- runif(1, 40, upper_bound)
  
}

dlog_prop_tau_PH <- function(tau_HP, d_obs_screen_i, d_obs_clinical_i, theta){
  
  age_screen_positive <- d_obs_screen_i %>% filter(screen_detected == 1) %>% pull(age_screen)
  upper_bound         <- min(age_screen_positive, d_obs_clinical_i$tau_PC_obs, na.rm = TRUE)
  dunif(tau_HP, 40, upper_bound, log = TRUE)
  
}

MH_tau_PH <- function(person_id, d_obs_screen_i, d_obs_clinical_i){
  
  # dat for person i
  d_obs_screen_i   <- d_obs_screen   %>% filter(person_id == i)
  d_obs_clinical_i <- d_obs_clinical %>% filter(person_id == i)
  d_process_i      <- d_process      %>% filter(person_id == i)
  
  # propose new tau_HP
  tau_PH_new <- rprop_tau_PH(d_obs_screen_i, d_obs_clinical_i)
  d_process_i_new            <- d_process_i
  d_process_i_new$tau_HP     <- tau_PH_new
  d_process_i_new$sojourn_H  <- tau_PH_new - d_process_i_new$age_start
  if(!d_process_i_new$indolent){
    d_process_i_new$sojourn_P <- d_process_i_new$tau_PC - d_process_i_new$tau_HP
  }
  
  # M-H acceptance ratio
  dlog_prop_new  <- dlog_prop_tau_PH(d_process_i_new$tau_HP, d_obs_screen_i, d_obs_clinical_i, theta)
  dlog_prop_curr <- dlog_prop_tau_PH(d_process_i$tau_HP    , d_obs_screen_i, d_obs_clinical_i, theta)
  dlog_lik_new   <- dlog_likelihood_i(d_process_i_new      , d_obs_screen_i, theta)
  dlog_lik_curr  <- dlog_likelihood_i(d_process_i          , d_obs_screen_i, theta)
  
  MH_logratio <- dlog_lik_new - dlog_lik_curr + dlog_prop_curr - dlog_prop_new
  
  if(runif(1) < MH_logratio){ # accept new tau_HP
    d_process[i, ] <- d_process_i_new
  }
  
  return(d_process)
}

#
# MCMC ####
MCMC <- function(){
  
  m <- 1e2
  RATE_H <- RATE_P <- BETA <- PSY <- numeric(m)
  
  for(n in 1 : m){
    
    # update parameters
    rate_H    <- gibbs_rate_H(d_process, theta, prior)
    rate_P    <- gibbs_rate_P(d_process, theta, prior)
    psy       <- gibbs_psy   (d_process, prior)
    rate_beta <- gibbs_beta  (d_process, d_obs_screen, prior)
    
    # update latent data
    for(i in 1 : nrow(d_process)){
      MH_tau_PH(i, )
      rprop_tau_PH(d_obs_screen_i, d_obs_clinical_i)
    }
    
    RATE_H[n] <- rate_H
    RATE_P[n] <- rate_P
    BETA  [n] <- psy
    PSY   [n] <- rate_beta
    
  }
  
  # output
  THETA <- list(RATE_H=RATE_H, RATE_P=RATE_P, BETA=BETA, PSY=PSY)
  out <- list(THETA)
  return(out)
  
}