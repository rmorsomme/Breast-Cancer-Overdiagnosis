i <- 2
indolent_i <- d_process$indolent[i]
tau_HP_i <- tau_HP[i]
d_obs_screen_i <- filter(d_obs_screen, person_id == i)
d_obs_censor_i <- filter(d_obs_censor, person_id == i)

censor_time_i <- d_obs_censor_i$censor_time
censor_type_i <- d_obs_censor_i$censor_type

taus <- seq(40, censor_time, 0.1)
m <- length(taus)
loglik <- sojourn_H <- sojourn_P <- screens <- numeric(m)


for(i in 1:m){ # pretty slow
  
  tau_HP <- taus[i]
  
  sojourn_H[i] <- dlog_sojourn_H   (tau_HP - 40         , theta)
  sojourn_P[i] <- dlog_sojourn_P   (censor_time - tau_HP, theta, indolent, d_obs_censor_i$censor_type)
  screens  [i] <- dlog_beta        (d_obs_screen_i, theta, tau_HP, censor_time)
  loglik   [i] <- dlog_likelihood_i(tau_HP, indolent, d_obs_screen_i, d_obs_censor_i, theta)
    
}

plot(taus, exp(sojourn_H))
plot(taus, exp(sojourn_P))
plot(taus, exp(sojourn_P + sojourn_H))
plot(taus, exp(screens  ))
plot(taus, exp(sojourn_H + screens  ))
plot(taus, exp(sojourn_P + screens  ))
plot(taus, exp(loglik))


####
ids <- which(d_obs_censor$censor_type=="censored")
i <- ids[1]
indolent_i <- d_process$indolent[i]
tau_HP_i <- tau_HP[i]
d_obs_screen_i <- filter(d_obs_screen, person_id == i)
d_obs_censor_i <- filter(d_obs_censor, person_id == i)

censor_time_i <- d_obs_censor_i$censor_time
censor_type_i <- d_obs_censor_i$censor_type

taus <- seq(40, censor_time, 0.1)
m <- length(taus)
loglik <- sojourn_H <- sojourn_P <- screens <- numeric(m)


for(i in 1:m){ # pretty slow
  
  tau_HP_i <- taus[i]
  
  sojourn_H[i] <- dlog_sojourn_H   (tau_HP_i, theta, censor_time_i)
  sojourn_P[i] <- dlog_sojourn_P   (tau_HP_i, theta, indolent_i, censor_type_i, censor_time_i)
  screens  [i] <- dlog_beta        (d_obs_screen_i, theta, tau_HP_i)
  loglik   [i] <- dlog_likelihood_i(tau_HP_i, indolent_i, d_obs_screen_i, d_obs_censor_i, theta)
  
}
