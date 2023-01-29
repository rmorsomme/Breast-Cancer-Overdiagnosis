
source("mcmc_setup.R")
load(file_id)

#
# Parameters ####

{
  THETA <- out$THETA
  
  RATE_H <- THETA$RATE_H#[(m/2):m]
  RATE_P <- THETA$RATE_P#[(m/2):m]
  PSI    <- THETA$PSI   #[(m/2):m]
  BETA   <- THETA$BETA  #[(m/2):m]
}

#
## Traceplot ####
plot(RATE_H, type = "l"); abline(h=theta$rate_H, col = "red")
plot(RATE_P, type = "l"); abline(h=theta$rate_P, col = "red")
plot(PSI   , type = "l"); abline(h=theta$psi   , col = "red")
plot(BETA  , type = "l"); abline(h=theta$beta  , col = "red")

#
## Histograms ####
x <- seq(1e-5, max(RATE_H, RATE_P, BETA, PSI), 0.0001)

{
  hist(RATE_H, freq=F, breaks=20)
  abline(v=theta$rate_H, col = "red")
  abline(v=quantile(RATE_H, c(0.025, 0.975)), col = "grey")
  lines(x, dgamma(x, prior$shape_H, prior$rate_H))
}

{
  hist(RATE_P, freq=F, breaks=20)
  abline(v=theta$rate_P, col = "red")
  abline(v=quantile(RATE_P, c(0.025, 0.975)), col = "grey")
  lines(x, dgamma(x, prior$shape_P, prior$rate_P))
}
{
  hist(PSI   , freq=F, breaks=20)
  abline(v=theta$psi   , col = "red")
  abline(v=quantile(PSI, c(0.025, 0.975)), col = "grey")
  lines(x, dbeta(x, prior$a_psi, prior$b_psi))
}
{
  hist(BETA  , freq=F, breaks=20)
  abline(v=theta$beta  , col = "red")
  abline(v=quantile(BETA, c(0.025, 0.975)), col = "grey")
  lines(x, dbeta(x, prior$a_beta, prior$b_beta))
}

#
## ACF ####
acf(RATE_H)
acf(RATE_P)
acf(PSI)
acf(BETA)

#
## Effective sample size ####
coda::effectiveSize(RATE_H)
coda::effectiveSize(RATE_P)
coda::effectiveSize(PSI)
coda::effectiveSize(BETA)

runtime <- out$runtime
coda::effectiveSize(RATE_H)/runtime
coda::effectiveSize(RATE_P)/runtime
coda::effectiveSize(PSI)   /runtime
coda::effectiveSize(BETA)  /runtime

#
## Correlation
THETA <- tibble(log(RATE_H), log(RATE_P), BETA, PSI)
round(cor(THETA), 2)
plot(THETA)


#
# Latent data ####
{
  TAU_HP       <- out$TAU_HP
  TAU_HP_inf   <- is.infinite(TAU_HP)
  ACCEPTED     <- out$ACCEPTED
  accept_rate  <- colMeans(ACCEPTED)
  tau_inf_rate <- colMeans(TAU_HP_inf)
}

#
## group ####

group_id <- c("clinical", "screen", "censored")[1]

{
  ids <- which(d_obs_censor$censor_type==group_id)
  accept_rate_group  <- accept_rate [ids]
  tau_inf_rate_group <- tau_inf_rate[ids]
  
  summary(accept_rate_group) %>% print
  order(accept_rate_group)[1:10] %>% print
}

if(group_id == "censored"){
  boxplot(tau_inf_rate_group)
  summary(tau_inf_rate_group) %>% print
  order(tau_inf_rate_group)[1:10] %>% print
}

if(group_id == "censored")  plot(d_obs_censor$censor_time[ids], tau_inf_rate[ids], xlab = "Censoring age (c)", ylab = "Probability of (c<tau)")
if(group_id == "censored")  plot(d_obs_censor$censor_time[ids], tau_inf_rate[ids], xlab = "Censoring age (c)", ylab = "Probability of (c<tau)", xlim=c(40,70))

boxplot(accept_rate_group)
if(group_id == "screen")  boxplot(accept_rate_group ~ d_obs_censor$censor_time[ids])
plot(d_obs_censor$censor_time[ids], accept_rate_group)

#
## individual ####
{
  i_lowest_accept  <- order(accept_rate_group)[1]
  i_highest_accept <- order(accept_rate_group)[length(accept_rate_group)]
  i <- ids[i_highest_accept]
  filter(d_obs_screen, person_id == i) %>% print
  d_process[i,] %>% print
  d_obs_censor[i,] %>% print
  
  mcmc_tau     <- TAU_HP[,i]
  mcmc_tau_inf <- TAU_HP_inf[,i]
  tau_true <- d_process[i, "tau_HP"]
}

accept_rate[i]

if(group_id == "censored")  plot(mcmc_tau_inf, type="l")

plot(mcmc_tau[is.finite(mcmc_tau)], type="l"); abline(h=tau_true, col = "red")

if(group_id == "censored"){
  acf(as.numeric(mcmc_tau_inf))
}else{
  acf(mcmc_tau)
}


if(group_id != "censored"){
  hist(mcmc_tau, breaks = 200, xlim=c(max(40, max(mcmc_tau) - 15), max(mcmc_tau)))
  abline(v=tau_true, col = "red")
  abline(v=quantile(mcmc_tau, c(0.025, 0.975)), col = "grey")
}

