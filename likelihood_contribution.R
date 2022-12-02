ids <- which(d_obs_censor$censor_type=="clinical")
i <- 866
indolent_i <- d_process$indolent[i]
tau_HP_i <- d_process$tau_HP[i]
d_obs_screen_i <- filter(d_obs_screen, person_id == i)
d_obs_censor_i <- filter(d_obs_censor, person_id == i)

censor_time_i <- d_obs_censor_i$censor_time
censor_type_i <- d_obs_censor_i$censor_type
endpoints_i <-  endpoints[[i]]

dt <- 1e-2
taus <- seq(max(40, censor_time_i - 15), censor_time_i, dt)
m <- length(taus)
loglik <- sojourn_H <- sojourn_P <- screens <- numeric(m)


for(t in 1:m){ # pretty slow
  
  tau_HP <- taus[t]
  
  sojourn_H[t] <- dlog_sojourn_H   (tau_HP, theta, censor_time_i)
  sojourn_P[t] <- dlog_sojourn_P   (tau_HP, theta, indolent_i, censor_type_i, censor_time_i)
  screens  [t] <- dlog_beta        (d_obs_screen_i, theta, tau_HP, censor_type_i)
  loglik   [t] <- dlog_likelihood_i(tau_HP, indolent_i, d_obs_screen_i, censor_type_i, censor_time_i, theta)
    
}

plot(taus, exp(sojourn_H))
plot(taus, exp(sojourn_P))
plot(taus, exp(sojourn_P + sojourn_H))
plot(taus, exp(screens  ))
plot(taus, exp(sojourn_H + screens  ))
plot(taus, exp(sojourn_P + screens  ))
plot(taus, exp(loglik))

