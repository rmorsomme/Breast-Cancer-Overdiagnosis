i <- 1
indolent <- d_process$indolent[i]
d_obs_screen_i <- filter(d_obs_screen, person_id == i)
d_obs_censor_i <- filter(d_obs_censor, person_id == i)

censor_time <- d_obs_censor_i$censor_time[i]
taus <- seq(55, censor_time, 0.1)
m <- length(taus)
loglik <- sojourn_H <- sojourn_P <- screens <- numeric(m)


for(i in 1:m){
  
  tau_HP <- taus[i]
  
  sojourn_H[i] <- dlog_sojourn_H   (tau_HP - 40         , theta)
  sojourn_P[i] <- dlog_sojourn_P   (censor_time - tau_HP, theta, indolent, d_obs_censor_i$censor_type)
  screens  [i] <- dlog_beta        (d_obs_screen_i, theta, tau_HP, censor_time)
  loglik   [i] <- dlog_likelihood_i(tau_HP, indolent, d_obs_screen_i, d_obs_censor_i, theta)
    
}

plot(taus, exp(sojourn_H))
plot(taus, exp(sojourn_P))
plot(taus, exp(screens  ))
plot(taus, exp(sojourn_P + sojourn_H))
plot(taus, exp(sojourn_H + screens  ))
plot(taus, exp(sojourn_P + screens  ))
plot(taus, exp(loglik))


m <- 1e5
tau_HP <- numeric(m)
for(i in 1:m){
  tau_HP[i] <- rprop_tau_HP(d_obs_screen_i, d_obs_censor_i, theta)
}
hist(tau_HP, breaks = 100)
mean(64<tau_HP)
mean(62<tau_HP&tau_HP<64)
round(prob, 3)
