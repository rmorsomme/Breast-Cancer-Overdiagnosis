{
    #
    # setup ####
    rm(list=ls())
    
    library(coda) # effectiveSize()
    library(tidyverse)
    
    setwd("C:/Users/18582/Desktop/Research/Marc/Breast Cancer Overdiagnosis")
    source("slurm/R/helpers.R")
    source("slurm/R/helpers-figures.R")
    
    # paths to files and figures
    path_mcmc  <- "slurm/output/MCMC/simulation"
    path_fig   <- "slurm/output/figures/simulation"
    
    # default configuration for figures
    scale = 6
    theme_set(theme_bw())
    #theme_update(text = element_text(size = 25))
    
    # initial values for MCMC
    theta_0_type = c("true", "low_density")[2]
    
    # Load data
    sim_id     <- paste0(
        "M=", 1e5,
        "-shape_H=", 2,
        "-shape_P=", 1,
        "-t0=", 30,
        "-theta_0=", theta_0_type,
        ".RDATA"
    )
    
    file_draws <- paste(path_mcmc, sim_id, sep = "/")
    
    load(file_draws)
    
    # Run time ####
    runtime <- out$runtime
    print(runtime / 3600) # hours
    
    # Acceptance rate ####
    print(mean(out$ACCEPT$ACCEPT_RATE_H))
    print(mean(out$ACCEPT$ACCEPT_RATE_P))
    print(mean(out$ACCEPT$ACCEPT_PSI   ))
}


#
# Data
d_obs_censor %>% count()


#
# Traceplot - transient


# Traceplot
burnin <- 500
out$THETA %>%
    as_tibble() %>%
    mutate(iteration = (1:nrow(.))*thin) %>%
    filter(iteration < burnin) %>%
    ggplot(aes(x = iteration, y = PSI)) +
    geom_line() +
    scale_x_continuous(breaks=seq(0,burnin,length.out=6)) +
    labs(x = "Iteration", y = expression(psi))
save_figures("traceplot_transient", path = path_fig, scale = scale)


#
# Traceplot - recurrent

# Traceplot
out$THETA %>%
    as_tibble() %>%
    mutate(
        iteration = (1:nrow(.))*thin,
        iteration_1000 = iteration / 1e3
    ) %>%
    filter(iteration > burnin) %>%
    ggplot(aes(x = iteration_1000, y = PSI)) +
    geom_line() +
    scale_x_continuous(breaks=seq(0,100,length.out=6)) +
    labs(x = "Iteration (in 1,000)", y = expression(psi))
save_figures("traceplot_recurrent", path = path_fig, scale = scale)


#
# ACF ###
my_acf <- function(draws){
    tmp <- acf(draws, plot = FALSE)
    tibble(acf = as.vector(tmp$acf), lag = as.vector(tmp$lag))
}

out$THETA %>%
    as_tibble() %>%
    pivot_longer(cols=RATE_H:PSI, names_to = "parameter", values_to = "draws") %>%
    mutate(
        parameter = case_when(
            parameter=="BETA"~"beta",
            parameter=="PSI"~"psi",
            parameter=="RATE_H"~"lambda[H]",
            parameter=="RATE_P"~"lambda[P]")
    ) %>%
    group_by(parameter) %>%
    summarize(test = list(my_acf(draws))) %>%
    unnest(test) %>%
    ggplot(aes(lag, acf)) + 
    geom_col() + 
    scale_y_continuous(breaks=seq(0,1,0.5)) +
    facet_wrap(~parameter, labeller = label_parsed) +
    labs(x = "Lags", y= "Auto-correlation\n function")
save_figures("acf", path = path_fig, scale = scale)

#
# summary ####
out$THETA %>%
    as_tibble() %>%
    pivot_longer(cols=RATE_H:PSI, names_to = "parameter", values_to = "draws") %>%
    group_by(parameter) %>%
    summarize(
        mean    = mean(draws),
        q_low   = quantile(draws, probs = 0.025),
        q_high  = quantile(draws, probs = 0.975),
        ESS     = effectiveSize(draws),
        ESS_sec = ESS / out$runtime
    )


#
# latent data

# find a screen-detected with multiple negative screens
d_obs_censor %>%
    filter(censor_type == "screen", censor_time-AFS>10)
id = 1708  
d_obs_screen_id = d_obs_screen %>% filter(person_id == id) %>% print()
id_col = which(d_obs_censor$person_id == id)

# indolent
indolent = out$INDOLENT
tibble(indolent = indolent[,id_col]) %>%
    mutate(
        iteration = (1:nrow(.))*thin,
        iteration_1000 = iteration / 1e3
    ) %>%
    ggplot(aes(iteration_1000, indolent)) +
    geom_line() +
    xlim(10,15) +
    labs(x = "Iteration (in 1,000)", y = "Latent\nindolent indicator")
save_figures("traceplot_indolent", path = path_fig, scale = scale)
id_col = which(d_obs_censor$person_id == id)

# tau
tau = out$age_at_tau_hp_hat
tau_tbl = tibble(tau = tau[,id_col]) %>%
    mutate(
        iteration = (1:nrow(.))*thin,
        iteration_1000 = iteration / 1e3
    )
tau_tbl %>%
    ggplot(aes(iteration_1000, tau)) +
    geom_line() +
    #xlim(10,15) +
    labs(x = "Iteration (in 1,000)", y = "Latent pre-clinial\ntransition age")
save_figures("traceplot_tau", path = path_fig, scale = scale)

tau_tbl %>%
    ggplot(aes(tau)) +
    geom_histogram(boundary=58, binwidth = 1/5) +
    scale_x_continuous(breaks=seq(0,100,by=1)) +
    geom_vline(xintercept = d_obs_screen_id$age_screen[10:15]) +
    labs(x = "Latent pre-clinial\ntransition age")
save_figures("histogram_tau", path = path_fig, scale = scale)


#
# Misc.
theta <- list(
    rate_H = 1/580, shape_H = 2,
    rate_P = 1/6  , shape_P = 1,
    beta   = 0.85, psi = 0.1
) %>% update_scales()

t0=30
age_screen = 40+5*(1:5)

for(censor_type in c("censored", "screen")) {
    proposal_vs_target(censor_type, age_screen, theta, t0)
}
