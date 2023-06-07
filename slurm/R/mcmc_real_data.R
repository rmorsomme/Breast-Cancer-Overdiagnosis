
rm(list = ls())

#
# Set up
#setwd("C:/Users/18582/Desktop/Research/Marc/Breast Cancer Overdiagnosis/slurm")
library(tidyverse)
library(tictoc)
library(parallel)
library(beepr)
source("R/helpers.R") # load helper functions
set.seed(0)

n_cpu <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))
if(is.na(n_cpu))  n_cpu <- 1


#
# Data ####
data_origin <- c("BCSC", "Swiss")[1]

if(data_origin == "BCSC"){
  
    load("data/processed/BCSC_40_to_85.RDATA")
    AFS_low <- 45
    AFS_upp <- 74
    n   <- min(1e5, nrow(d_obs_censor %>% filter(AFS %>% between(AFS_low, AFS_upp))))
  
    d_obs_censor <- d_obs_censor %>%
        filter(AFS %>% between(AFS_low, AFS_upp)) %>%
        slice_sample(n = n) %>%
        arrange(person_id)
  
    d_obs_screen <- d_obs_screen %>%
        filter(person_id %in% d_obs_censor$person_id)
    
}else if(data_origin == "Swiss"){
  
    load("data/processed/swiss.RDATA")
    AFS_low <- NA
    AFS_upp <- NA
  
}

  


#
# Model setup ####
t0      <- 30 # for Swiss data, need to check that there is no screen before age t0
shape_H <- 2  # should be 2 (linear hazard rate) or larger
shape_P <- 1


#
# MCMC setup ####
M       <- 1e4 # number of MCMC iterations
thin    <- round(max(M/1e3, 1))

epsilon_rate_H <- 2e-5 # step size of RW MH
epsilon_rate_P <- 0.1
epsilon_psi    <- 0.075

mean_H_0 <- 85
rate_H_0 <- (mean_H_0/gamma(1+1/shape_H))^(-shape_H) # this gives a Weibull with mean mean_H_0
theta_0  <- list( # initial values for theta
    rate_H = rate_H_0, shape_H = shape_H,
    rate_P = 0.2  , shape_P = shape_P,
    beta   = 0.8  , psi = 0.1
    ) %>%
    update_scales()


a <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")) - 1
if(is.na(a))  a <- 0

#precision_mean_P <- c(3, 15, 100, 1000)[1 + a %% 4]
#mean_mean_P      <- c(3, 5, 7)[1 + a %/% 4]
precision_mean_P <- 0
mean_mean_P      <- 0

prior <- list(
    #rate_H = 0.01, shape_H = 1, # gamma(shape_H, rate_H) prior on the Weibull rate for H
    #rate_P = mean_mean_P * (precision_mean_P - 1), shape_P = precision_mean_P, # gamma(shape_P, rate_P) prior on the Weibull (exponential) rate for P
    rate_H = 0.01, shape_H = 1,
    rate_P = 0.01, shape_P = 1,
    a_psi  = 1   , b_psi   = 1,   # beta(a_psi , b_psi ) prior on psi
    a_beta = 38.5, b_beta  = 5.8  # beta(a_beta, b_beta) prior on beta
)


#
# MCMC ####
out <- MCMC(
    d_obs_screen, d_obs_censor,
    theta_0, prior,
    epsilon_rate_H, epsilon_rate_P, epsilon_psi, 
    t0, M, thin, n_cpu,
    verbose = TRUE
    )


#
# Output ####

path_mcmc  <- paste0("output/MCMC/", data_origin)
sim_id     <- paste0(
    "M=", M,
    "-AFS_low=", AFS_low,
    "-AFS_upp=", AFS_upp,
    "-shape_H=", shape_H,
    "-shape_P=", shape_P,
    "-t0=", t0,
    "-mean_mean_P=", mean_mean_P,
    "-precision_mean_P=", precision_mean_P,
    ".RDATA"
)
file_draws <- paste(path_mcmc, sim_id, sep = "/")

theta <- NULL

save(
    out,
    M, thin, n_cpu,
    theta, 
    prior, theta_0, epsilon_rate_H, epsilon_rate_P, epsilon_psi,
    d_obs_screen, d_obs_censor,
    file = file_draws
)
