loglik <- sojourn_H <- sojourn_P <- screens <- c()

taus <- seq(55, 66, 0.1)

for(tau in taus){
  
  tau_PH_new <- tau
  d_process_i_new            <- d_process_i
  d_process_i_new$tau_HP     <- tau_PH_new
  d_process_i_new$sojourn_H  <- tau_PH_new - d_process_i_new$age_start
  if(!d_process_i_new$indolent){
    d_process_i_new$sojourn_P <- d_process_i_new$tau_PC - d_process_i_new$tau_HP
  }
  
  loglik    <- c(loglik   , dlog_likelihood_i(d_process_i_new , d_obs_screen_i, theta))
  sojourn_H <- c(sojourn_H, dlog_sojourn_H(d_process_i_new$sojourn_H, theta))
  sojourn_P <- c(sojourn_P, dlog_sojourn_P(d_process_i_new$sojourn_P, theta))
  screens   <- c(screens  , dlog_beta(d_obs_screen_i, theta, d_process_i_new$tau_HP, d_process_i_new$tau_PC))
    
  
}

plot(taus, exp(loglik))
plot(taus, exp(sojourn_H))
plot(taus, exp(sojourn_P))
plot(taus, exp(screens))


