
rm(list = ls())

#
# Set up
library(tidyverse)
library(tictoc)
source("helpers.R") # load helper functions
set.seed(0)

theta <- list(
    rate_H = 0.0005, shape_H = 2,
    rate_P = 0.0015, shape_P = 4,
    beta   = 0.814 , psi     = 0.75 # psi = 0.045
)
theta <- update_scales(theta)

#
# Simulate data
n <- 1e3  # sample size
source("simulator.R")

#
# MCMC setup
m <- 2e3 # number of MCMC iterations

factor_0 <- 1
theta_0  <- list(
  rate_H = theta$rate_H*factor_0, shape_H = theta$shape_H,
  rate_P = theta$rate_P*factor_0, shape_P = theta$shape_P,
  beta   = min(1-0.01, theta$beta*factor_0), psi = min(1-0.01, theta$psi*factor_0)
)
theta_0  <- update_scales(theta_0)

prior <- list(
  rate_H = 1, shape_H = 1, # gamma(shape_H, rate_H) prior on Weibull rate for H
  rate_P = 1.333e3, shape_P = 2, # gamma(shape_P, rate_P) prior on Weibull rate for P
  a_psi  = 1, b_psi   = 1, # beta(a_psi, b_psi)     prior on psi
  a_beta = 8, b_beta  = 2  # beta(a_beta, b_beta)   prior on beta
)

# path to save mcmc draws
path_output <- "output"
path_mcmc   <- paste0(path_output, "/MCMC")
file_id     <- paste0(path_mcmc, "/MCMC-n=", n, "-m=", m, "-factor0=", factor_0, ".RDATA")
