#
# Set up
#setwd("C:/Users/18582/Desktop/Research/Marc/Breast Cancer Overdiagnosis/slurm")
rm(list = ls())
set.seed(0)

library(tidyverse)
library(tictoc)
library(beepr)
library(Rcpp)

sourceCpp("R/helpers_in_Rcpp.cpp")
source("R/helpers_for_calling_Rcpp.R")

#
# load artificial data ####
path_data  <- "output/MCMC/simulation/data"
id_data     <- paste0(
    "n=", 37e3,
    "-shape_H=", 2,
    "-shape_P=", 1,
    "-t0=", 30,
    ".RDATA"
)

file_data <- paste(path_data, id_data, sep = "/")
load(file_data)


#
# MCMC setup ####
M       <- 1e5 # total number of MCMC iterations
M_thin  <- 1e4 # number of MCMC iterations kept
thin    <- round(max(M/M_thin, 1))

# tuning parameter
epsilon_rate_H <- 9e-6
epsilon_rate_P <- 0.05
epsilon_psi    <- 0.1

# prior distribution
prior <- list(
    rate_H = 0.01, shape_H = 1, # gamma(shape_H, rate_H) prior on Weibull rate for H
    rate_P = 0.01, shape_P = 1, # gamma(shape_H, rate_H) prior on Weibull rate for P
    a_psi  = 1/2 , b_psi   = 1/2, # beta(a_psi , b_psi )   prior on psi
    a_beta = 38.5, b_beta  = 5.8  # beta(a_beta, b_beta)   prior on beta
)

# initial values, either the true values or values in a low density region
theta_0_type = c("true", "low_density")[1]
theta_0 = if(theta_0_type == "true") {
    theta
} else if(theta_0_type == "low_density") {
    list(
        rate_H = theta$rate_H, shape_H = theta$shape_H,
        rate_P = theta$rate_P, shape_P = theta$shape_P,
        beta   = theta$beta, psi = 0.8
    ) %>% update_scales()
}

#
# MCMC ####
out <- MCMC_cpp(
    d_obs_screen, d_obs_censor,
    theta_0, prior, 
    epsilon_rate_H, epsilon_rate_P, epsilon_psi, 
    t0, M, thin
)


#
# Save output ####
path_mcmc  <- "output/MCMC/simulation"
path_fig   <- "output/figures/simulation"
sim_id     <- paste0(
    "M=", M,
    "-shape_H=", shape_H,
    "-shape_P=", shape_P,
    "-t0=", t0,
    "-theta_0=", theta_0_type,
    ".RDATA"
)

# Save parameters only
file_draws <- paste(path_mcmc, sim_id, sep = "/")
save(
    out,
    M, thin,
    theta,
    prior, theta_0, epsilon_rate_H, epsilon_rate_P, epsilon_psi,
    d_obs_screen, d_obs_censor,
    file = file_draws
)


#
# Acceptance rate ####
print(mean(out$ACCEPT$ACCEPT_RATE_H))
print(mean(out$ACCEPT$ACCEPT_RATE_P))
print(mean(out$ACCEPT$ACCEPT_PSI   ))

beep()