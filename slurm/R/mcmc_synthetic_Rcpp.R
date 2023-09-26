
#
# Set up
setwd("C:/Users/18582/Desktop/Research/Marc/Breast Cancer Overdiagnosis/slurm")
rm(list = ls())

library(tidyverse)
library(tictoc)
library(parallel)
library(beepr)
library(Rcpp)

source("R/helpers.R")
sourceCpp("R/helpers_in_Rcpp.cpp")
source("R/helpers_for_calling_Rcpp.R")

set.seed(1)

n_cpu <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))
if(is.na(n_cpu))  n_cpu <- 1

#
# Model setup ####
t0 <- 30
shape_H <- 2 # should be 2 (linear) or larger
shape_P <- 1


#
# Parameters ####
theta <- list(
    rate_H = 1/500, shape_H = shape_H,
    rate_P = 1/5  , shape_P = shape_P,
    beta   = 38.5/(38.5+5.8), psi = 0.1
) %>% update_scales()


#
## Data ####
n  <- 1e4 #37e3  # sample size
source("R/simulator.R")


#
# MCMC setup ####
M       <- 1e4 # number of MCMC iterations
thin    <- round(max(M/1e3, 1))

epsilon_rate_H <- 0.001 # tuning parameter
epsilon_rate_P <- 0.05
epsilon_psi    <- 0.1

factor_0 <- 1
theta_0  <- list(
    rate_H = theta$rate_H*factor_0, shape_H = theta$shape_H,
    rate_P = theta$rate_P*factor_0, shape_P = theta$shape_P,
    beta   = min(1-0.01, theta$beta*factor_0), psi = min(1-0.01, theta$psi*factor_0)
) %>% update_scales()


precision_mean_P <- 1
mean_mean_p      <- 5
prior <- list(
    rate_H = 2e1 , shape_H = 2, # gamma(shape_H, rate_H) prior on Weibull rate for H
    rate_P = mean_mean_p * precision_mean_P , shape_P = precision_mean_P + 1, # gamma(shape_P, rate_P) prior on Weibull (exponential) rate for P
    a_psi  = 1   , b_psi   = 1,   # beta(a_psi , b_psi )   prior on psi
    a_beta = 38.5, b_beta  = 5.8  # beta(a_beta, b_beta)   prior on beta
)


#
# MCMC ####
out_cpp <- MCMC_cpp(
    d_obs_screen, d_obs_censor,
    theta_0, prior, 
    epsilon_rate_H, epsilon_rate_P, epsilon_psi, 
    t0, M, thin
)
out_r <- MCMC(
    d_obs_screen, d_obs_censor,
    theta_0, prior, 
    epsilon_rate_H, epsilon_rate_P, epsilon_psi, 
    t0, M, thin,
    n_cpu, verbose = TRUE
)

out_r$runtime / out_cpp$runtime

# #
# # Output ####
# path_mcmc  <- "output/MCMC/simulation"
# sim_id     <- paste0(
#     "n=", n, 
#     "-M=", M,
#     "-shape_H=", shape_H, 
#     "-shape_P=", shape_P,
#     ".RDATA"
# )
# file_draws <- paste(path_mcmc, sim_id, sep = "/")

beep()
