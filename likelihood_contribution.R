ids <- which(d_obs_censor$censor_type=="screen")
#ids <- which(d_obs_censor$censor_type=="clinical")
i <- ids[1]
indolent_i <- d_process$indolent[i]
indolent_i <- 1
tau_HP_i <- d_process$tau_HP[i]
d_obs_screen_i <- filter(d_obs_screen, person_id == i)
d_obs_censor_i <- filter(d_obs_censor, person_id == i)

censor_time_i       <- d_obs_censor_i$censor_time
censor_type_i       <- d_obs_censor_i$censor_type
endpoints_i         <- endpoints[[i]]
age_screen_i        <- age_screen[[i]]
n_screen_positive_i <- n_screen_positive[[i]]

dt <- 1e-2
taus <- seq(max(40 + 0.01, censor_time_i - 30), censor_time_i - 0.01, dt)
m <- length(taus)
loglik <- sojourn_H <- sojourn_P <- screens <- numeric(m)


for(t in 1:m){ # pretty slow
  
  tau_HP <- taus[t]
  
  sojourn_H[t] <- dlog_sojourn_H   (tau_HP, theta, censor_time_i)
  sojourn_P[t] <- dlog_sojourn_P   (tau_HP, theta, indolent_i, censor_type_i, censor_time_i)
  screens  [t] <- dlog_beta        (age_screen_i, theta, tau_HP, n_screen_positive_i)
  indolent [t] <- dlog_psi         (indolent_i, theta)
  loglik   [t] <- dlog_likelihood_i(
    tau_HP, indolent_i, censor_type_i, censor_time_i, age_screen_i, n_screen_positive_i, theta
    )
    
}

plot(taus, exp(sojourn_H), type = "l", ylim=c(0,max(exp(sojourn_H))))
plot(taus, exp(sojourn_P), type = "l")
plot(taus, exp(sojourn_P + sojourn_H), type = "l")
plot(taus, exp(screens  ), type = "l")
plot(taus, exp(sojourn_H + screens  ), type = "l")
plot(taus, exp(sojourn_P + screens  ), type = "l")
plot(taus, exp(loglik), type = "l")


normalize_dens <- function(f){
  f/sum()
  #n <- length(f)
  #f / sum((f[1:(n-1)] + f[2:n])/2 * diff(x))
}

# proposal versus full conditional
tibble(tau_HP = taus) %>%
  mutate(
    posterior_unnormalize = tau_HP %>% 
      map_dbl(dlog_likelihood_i, indolent_i, censor_type_i, censor_time_i, age_screen_i, n_screen_positive_i, theta) %>% 
      exp(),
    posterior = posterior_unnormalize / sum(posterior_unnormalize)
    )

,
    proposal  = tau_HP %>% 
      map_dbl(dlog_prop_tau_HP, censor_type_i, censor_time_i, endpoints_i, prob, theta) %>%
      exp() %>% 
      normalize_dens(tau_HP, taus)
    ) %>%
  pivot_longer(cols = posterior:proposal, names_to = "Distribution", values_to = "Density") %>%
  ggplot(aes(tau_HP, Density, col = Distribution)) +
  geom_line() +
  labs(x=expression(tau[i]^HP), y= "density") #expression( pi(tau[866]^HP  ~ "|" ~  theta) )

