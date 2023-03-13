
rm(list = ls())

#
# Set up
library(tidyverse)
library(tictoc)
library(parallel)
library(beepr)
source("helpers.R") # load helper functions
set.seed(0)

n_cpu <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))
if(is.na(n_cpu)) n_cpu <- 1

theta <- list(
    rate_H = 0.0005, shape_H = 2,
    rate_P = 0.0015, shape_P = 4,
    beta   = 0.814 , psi     = 0.1 # psi = 0.045
)
theta <- update_scales(theta)

#
# Simulate data
n <- 1e2  # sample size
source("simulator.R")

#
# MCMC setup
M    <- 1e3 # number of MCMC iterations
thin <- round(max(M/1e3, 1))

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
  a_psi  = 1, b_psi  = 1, # beta(a_psi, b_psi)     prior on psi
  a_beta = theta$beta*1e1 , b_beta = (1-theta$beta)*1e1  # beta(a_beta, b_beta)   prior on beta
)

epsilon <- 0.5
# path to save mcmc draws
path_mcmc   <- "output/MCMC"
file_id     <- paste0(path_mcmc, "/MCMC-n=", n, "-M=", M, "-factor0=", factor_0, "-n_cpu=", n_cpu, ".RDATA")
