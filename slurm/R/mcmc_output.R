

{
    rm(list = ls())
    #setwd("C:/Users/18582/Desktop/Research/Marc/Breast Cancer Overdiagnosis/slurm")
    library(tidyverse)
    #source("R/helpers.R") # load helper functions
    #load("data/processed/BCSC_40_to_85.RDATA")
    
    M <- 1e4
    
    data_origin <- c("BCSC", "Swiss", "simulation")[1]
    path_mcmc  <- paste0("output/MCMC/"   , data_origin)
    path_fig   <- paste0("output/figures/", data_origin)
    shape_H <- 2 # 1, 1.1, 1.5, 2 ,2.5, 3, 3.5
    shape_P <- 1
    #AFS_low <- 45
    #AFS_upp <- 74
    t0  <- 30
    precision_mean_P <- 0
    mean_mean_P <- 0
    
    sim_id     <- paste0(
        "M=", M,
        "-AFS_low=", AFS_low,
        "-AFS_upp=", AFS_upp,
        "-shape_H=", shape_H,
        "-shape_P=", shape_P,
        "-t0=", t0,
        #"-mean_mean_P=", mean_mean_P,
        #"-precision_mean_P=", precision_mean_P,
        ".RDATA"
    )
    file_draws <- paste(path_mcmc, sim_id, sep = "/")
    file_fig   <- paste(path_fig , sim_id, sep = "/")
    load(file_draws)
    
    theta$mean_H <- theta$rate_H^(-1/shape_H)*gamma(1+1/shape_H)
    theta$mean_P <- theta$rate_P^(-1/shape_P)*gamma(1+1/shape_P)
   
}

{
#
# Run time ####
(runtime <- out$runtime) / 60 # minutes
(runtime <- out$runtime) / 3600 # hours


#
# Acceptnace rate ####
mean(out$ACCEPT$ACCEPT_RATE_H)
mean(out$ACCEPT$ACCEPT_RATE_P)
mean(out$ACCEPT$ACCEPT_PSI   )



#
# Parameters ####

THETA <- out$THETA %>%
    as_tibble() %>%
    mutate(
        MEAN_H = RATE_H^(-1/theta_0$shape_H)*gamma(1+1/theta_0$shape_H),
        MEAN_P = RATE_P^(-1/theta_0$shape_P)*gamma(1+1/theta_0$shape_P),
        iteration = 1:nrow(.)
    ) %>%
    #filter(iteration > 100) %>% # remove burnin
    pivot_longer(cols = -iteration, names_to = "parameter", values_to = "draws")


{
    g <- ggplot(THETA, aes(iteration, draws)) +
        geom_line()+
        facet_wrap(~parameter, scales = "free") +
        labs(title=paste0("AFS between ", AFS_low, " and ", AFS_upp))
    print(g)
    ggsave(
        paste(sim_id, "traceplot.jpg", sep = "_"),
        path = path_fig, width = 1.61803, height = 1, scale = 5
    )
}
}
THETA_summary <- THETA %>%
    group_by(parameter) %>%
    summarize(
        mean    = mean(draws),
        q_low   = quantile(draws, probs = 0.025),
        q_high  = quantile(draws, probs = 0.975),
        ESS     = coda::effectiveSize(draws),
        ESS_sec = ESS / out$runtime
    ) %>%
    print()

my_acf <- function(draws){
    tmp <- acf(draws, plot = FALSE)
    tibble(
        acf = as.vector(tmp$acf),
        lag = as.vector(tmp$lag)
    )
}

THETA_acf <- THETA %>%
    group_by(parameter) %>%
    summarize(test = list(my_acf(draws))) %>%
    unnest(test)

THETA_true <- as_tibble(theta)

{
    g <- ggplot(THETA, aes(draws)) +
        geom_histogram(aes(y = after_stat(density))) +
        facet_wrap(~parameter, scales = "free") +
        geom_vline(data = THETA_summary, mapping = aes(xintercept=q_low)) +
        geom_vline(data = THETA_summary, mapping = aes(xintercept=q_high))
    print(g)
    ggsave(
        paste(sim_id, "histogram_raw.jpg", sep = "_"),
        path = path_fig, width = 1.61803, height = 1, scale = 5
    )
}

{
    g <- ggplot(THETA_acf, aes(lag, acf)) + 
        geom_col() + 
        facet_wrap(~parameter)
    print(g)
    ggsave(
        paste(sim_id, "acf.jpg", sep = "_"),
        path = path_fig, width = 1.61803, height = 1, scale = 5
    )
}


#
# Hazard rate
{
    age <- seq(t0, 75, by = 0.01)
    shape_H <- theta_0$shape_H
    rate_H_upp  <- THETA_summary$q_high[5]
    rate_H_low  <- THETA_summary$q_low [5]
    rate_H_mean <- THETA_summary$mean  [5]
    hazard_upp  <- shape_H * rate_H_upp  * (age-t0)^(shape_H-1)
    hazard_mean <- shape_H * rate_H_low  * (age-t0)^(shape_H-1)
    hazard_low  <- shape_H * rate_H_mean * (age-t0)^(shape_H-1)
    
    plot(age, hazard_upp, type = "l", main = "Hazard rate (posterior mean and 95% credible band)")
    lines(age, hazard_low)
    lines(age, hazard_mean)
    
    
    # hazard at ages 47.5, 60 and 70
    print(signif(shape_H * c(rate_H_low, rate_H_mean, rate_H_upp) * (47.5 - t0)^(shape_H-1)), 3)
    print(signif(shape_H * c(rate_H_low, rate_H_mean, rate_H_upp) * (60   - t0)^(shape_H-1)), 3)
    print(signif(shape_H * c(rate_H_low, rate_H_mean, rate_H_upp) * (70   - t0)^(shape_H-1)), 3)
}

#
## Correlation
out$THETA %>%
    as_tibble() %>%
    mutate(
        MEAN_H = RATE_H^(-1/theta_0$shape_H)*gamma(1+1/theta_0$shape_H),
        MEAN_P = RATE_P^(-1/theta_0$shape_P)*gamma(1+1/theta_0$shape_P),
        iteration = 1:nrow(.)
    ) %>%
    filter(iteration > 250) %>%
    cor() %>%
    round(2)

out$THETA %>%
    as_tibble() %>%
    mutate(
        MEAN_H = RATE_H^(-1/theta_0$shape_H)*gamma(1+1/theta_0$shape_H),
        MEAN_P = RATE_P^(-1/theta_0$shape_P)*gamma(1+1/theta_0$shape_P),
        iteration = 1:nrow(.)
    ) %>%
    #filter(iteration > 500) %>%
    select(-iteration) %>%
    GGally::ggpairs()



#
# Latent data ####
{
    TAU_HP             <- out$age_at_tau_hp_hat
    TAU_HP_inf         <- is.infinite(TAU_HP)
    ACCEPT_LATENT      <- out$ACCEPT$ACCEPT_LATENT
    ACCEPT_PSI         <- out$ACCEPT$ACCEPT_PSI
    tau_inf_rate       <- colMeans(TAU_HP_inf)
    accept_latent_rate <- colMeans(ACCEPT_LATENT)
    
    mean(ACCEPT_PSI) %>% print()
    mean(accept_latent_rate) %>% print()
}
#
## group ####

group_id <- c("clinical", "screen", "censored")[1]

{
    id_group <- which(d_obs_censor$censor_type==group_id)
    accept_rate_group  <- accept_latent_rate[id_group]
    tau_inf_rate_group <- tau_inf_rate[id_group]
    
    summary(accept_rate_group) %>% print()
}

if(group_id == "censored"){
    boxplot(tau_inf_rate_group)
    summary(tau_inf_rate_group) %>% print
    order(tau_inf_rate_group)[1:10] %>% print
}

if(group_id == "censored")  plot(d_obs_censor$censor_time[id_group], tau_inf_rate[id_group], xlab = "Censoring age (c)", ylab = "Probability of (c<tau)")
if(group_id == "censored")  plot(d_obs_censor$censor_time[id_group], tau_inf_rate[id_group], xlab = "Censoring age (c)", ylab = "Probability of (c<tau)", xlim=c(40,50))

boxplot(accept_rate_group)
if(group_id == "censored")  plot(tau_inf_rate_group, accept_rate_group)
if(group_id == "screen")  boxplot(accept_rate_group ~ d_obs_censor$censor_time[id_group])
plot(d_obs_censor$censor_time[id_group], accept_rate_group)
plot(d_obs_censor$censor_time[id_group], accept_rate_group, xlim = c(42, 52))

#
## individual ####
{
    i_lowest_accept  <- order(accept_rate_group)[1]
    i_highest_accept <- order(accept_rate_group)[length(accept_rate_group)]
    i <- id_group[i_lowest_accept]
    accept_latent_rate[i]
    
    d_obs_censor[i,] %>% print()
    
    hist(TAU_HP[,i])
    plot(TAU_HP[,i])
}


if(group_id == "censored")  plot(TAU_HP_inf[,i], type="l")

#plot(mcmc_tau[is.finite(mcmc_tau)], type="l"); abline(h=tau_true, col = "red")

if(group_id == "censored"){
    acf(as.numeric(TAU_HP_inf[,i]))
}else{
    acf(TAU_HP[,i])
}


if(group_id != "censored"){
    hist(mcmc_tau, breaks = 200, xlim=c(max(40, max(mcmc_tau) - 15), max(mcmc_tau)))
    abline(v=tau_true, col = "red")
    abline(v=quantile(mcmc_tau, c(0.025, 0.975)), col = "grey")
}
