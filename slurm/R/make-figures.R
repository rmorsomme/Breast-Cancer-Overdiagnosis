
#
# Misc. ####
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


#
# Simulations ####

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
d_obs_censor %>% count(censor_type)
d_obs_screen %>% count()


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
# latent data ####

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
# BCSC ####
{
    source("slurm/R/helpers.R")
    source("slurm/R/helpers-figures.R")
    scale = 6
    
    # paths to files and figures
    path_mcmc  <- "slurm/output/MCMC/BCSC"
    path_fig   <- "slurm/output/figures/BCSC"
    
    # Load data
    shape_P = 1
    shape_H = 2
    t0 = 30
    sim_id     <- paste0(
        "M=", 1e5,
        "-AFS_low=", 50,
        "-AFS_upp=", 74,
        "-shape_H=", shape_H,
        "-shape_P=", shape_P,
        "-t0=", t0,
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

# data
d_obs_censor %>% count()
d_obs_censor %>% count(censor_type)
d_obs_screen %>% count()


# marginal posterior: psi, mean_P
tbl_theta = out$THETA %>%
    as_tibble() %>%
    mutate(
        MEAN_P = RATE_P^(-1/shape_P)*gamma(1+1/shape_P),
        MEAN_H = RATE_H^(-1/shape_H)*gamma(1+1/shape_H)
        )

tbl_theta %>%
    ggplot(aes(x = PSI)) +
    geom_density() +
    geom_vline(xintercept = mean(tbl_theta$PSI), linetype="dashed", colour = "maroon") +
    geom_vline(xintercept = quantile(tbl_theta$PSI, c(0.025, 0.975)), linetype="dotted", colour = "maroon") +
    labs(x = expression(psi), y = "Density")
save_figures("histogram_psi", path = path_fig, scale = scale)

tbl_theta %>%
    as_tibble() %>%
    ggplot(aes(x = MEAN_P)) +
    geom_density() +
    geom_vline(xintercept = mean(tbl_theta$MEAN_P), linetype="dashed", colour = "maroon") +
    geom_vline(xintercept = quantile(tbl_theta$MEAN_P, c(0.025, 0.975)), linetype="dotted", colour = "maroon") +
    labs(x = "Expected pre-clinical\nsojourn time", y = "Density")
save_figures("histogram_mean_P", path = path_fig, scale = scale)

# hazard
rate_H_quantile = quantile(tbl_theta$RATE_H, c(0.025, 0.5, 0.975))
tibble(ages = seq(t0, 90, by = 0.1)) %>%
    mutate(
        h_low = shape_H * rate_H_quantile[1] * (ages-t0)^(shape_H-1),
        h_med = shape_H * rate_H_quantile[2] * (ages-t0)^(shape_H-1),
        h_upp = shape_H * rate_H_quantile[3] * (ages-t0)^(shape_H-1),
        ) %>%
    #pivot_longer(cols=h_low:h_upp, names_to = "hazard_type", values_to = "hazard") %>%
    ggplot(aes(x=ages)) +
    geom_line(aes(y=h_med)) +
    geom_line(aes(y=h_low), linetype = "dashed") +
    geom_line(aes(y=h_upp), linetype = "dashed") +
    labs(x="Age", y = "Hazard rate for\npre-clinical cancer")
save_figures("hazard", path = path_fig, scale = scale)

# hazard at ages 47.5, 60 and 70
print(signif(shape_H * rate_H_quantile * (47.5 - t0)^(shape_H-1)), 2)
print(signif(shape_H * rate_H_quantile * (60   - t0)^(shape_H-1)), 2)
print(signif(shape_H * rate_H_quantile * (70   - t0)^(shape_H-1)), 2)

# posterior summary
tbl_theta %>%
    pivot_longer(cols = RATE_H:MEAN_H, names_to = "parameter", values_to = "draws") %>%
    group_by(parameter) %>%
    summarize(
        mean    = mean(draws),
        q_low   = quantile(draws, probs = 0.025),
        q_high  = quantile(draws, probs = 0.975),
        ESS     = coda::effectiveSize(draws),
        ESS_sec = ESS / out$runtime
    )
