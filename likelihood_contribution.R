
group_id <- c("clinical", "screen", "censored")[1]
ids <- which(d_obs_censor$censor_type==group_id)
i <- ids[1]
indolent_i <- d_process$indolent[i]
#indolent_i <- 1
tau_HP_i <- d_process$tau_HP[i]
d_obs_screen_i <- filter(d_obs_screen, person_id == i)
d_obs_censor_i <- filter(d_obs_censor, person_id == i)

censor_time_i       <- d_obs_censor_i$censor_time
censor_type_i       <- d_obs_censor_i$censor_type
endpoints_i         <- endpoints[[i]]
age_screen_i        <- age_screen[[i]]
n_screen_positive_i <- n_screen_positive[[i]]

dt <- 1e-2
taus <- seq(max(40 + 0.01, censor_time_i - 10), censor_time_i - 0.01, dt)
m <- length(taus)
loglik <- sojourn_H <- sojourn_P <- screens <- numeric(m)


for(t in 1:m){
  
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
plot(taus, exp(sojourn_P), type = "l", ylim=c(0,max(exp(sojourn_P))))
plot(taus, exp(sojourn_P + sojourn_H), type = "l")
plot(taus, exp(screens  ), type = "l")
plot(taus, exp(sojourn_H + screens  ), type = "l")
plot(taus, exp(sojourn_P + screens  ), type = "l")
plot(taus, exp(loglik), type = "l")

tibble(tau = taus, loglik=loglik) %>%
    ggplot(aes(tau, loglik)) + 
    geom_line() +
    xlab(expression(t[i]))

ggsave("output/images/loglik.jpeg", height = 3.5, width = 3.5*1.615)



plot(taus, exp(loglik), type = "l")

prob_tau_i <- compute_prob_tau_i(indolent_i, censor_type_i, censor_time_i, endpoints_i, theta)
plot(prob_tau_i)
prob_tau_i <- prob_tau_i + c(0, -0.05, 0.05)

# proposal versus full conditional
tibble(tau_HP = taus) %>%
  mutate(
    posterior = tau_HP %>% 
      map_dbl(dlog_likelihood_i, indolent_i, censor_type_i, censor_time_i, age_screen_i, n_screen_positive_i, theta) %>% 
      exp(),
    posterior = posterior / sum(posterior),
    proposal  = tau_HP %>% 
      map_dbl(dlog_prop_tau_HP_i, indolent_i, censor_type_i, censor_time_i, endpoints_i, prob_tau_i, theta) %>%
      exp(),
    proposal = proposal / sum(proposal)
    ) %>%
  pivot_longer(cols = posterior:proposal, names_to = "Distribution", values_to = "Density") %>%
  ggplot(aes(tau_HP, Density, col = Distribution)) +
  geom_line() +
  labs(x=expression(tau[i]^HP), y= "density") #expression( pi(tau[866]^HP  ~ "|" ~  theta) )
