
rm(list = ls())

#
# Set up

n_cpu <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))
if(is.na(n_cpu)){
  n_cpu <- 1
  setwd("C:/Users/18582/Desktop/Research/Marc/Breast Cancer Overdiagnosis/slurm")
}

library(tidyverse)
library(tictoc)
library(parallel)
library(beepr)
library(Rcpp)

#source("R/helpers.R")
sourceCpp("R/helpers_in_Rcpp.cpp")
source("R/helpers_for_calling_Rcpp.R")

set.seed(1)


#
# Data ####
data_origin <- c("BCSC", "Swiss")[1]

if(data_origin == "BCSC") {
  
    load("data/processed/BCSC_40_to_85.RDATA")
    AFS_low <- 40
    AFS_upp <- 49
    #n   <- min(1e5, nrow(d_obs_censor %>% filter(AFS %>% between(AFS_low, AFS_upp))))
    
    #plot(table(d_obs_censor$AFS))
  
    d_obs_censor <- d_obs_censor %>%
        filter(AFS %>% between(AFS_low, AFS_upp)) %>%
        slice_sample(n = 1e4) %>%
        arrange(person_id)
  
    d_obs_screen <- d_obs_screen %>%
        filter(person_id %in% d_obs_censor$person_id)
} else if(data_origin == "Swiss") {
  
    load("data/processed/swiss.RDATA")
    AFS_low <- NA
    AFS_upp <- NA
}

  

#
# Model setup ####
t0      <- 30
shape_H <- 2
shape_P <- 1


#
# MCMC setup ####
M       <- 1e4 # number of MCMC iterations
thin    <- round(max(M/1e3, 1))

# step size of RW
epsilon_rate_H <- 5e-6
epsilon_rate_P <- 0.01
epsilon_psi    <- 0.1

# a <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")) - 1
# if(is.na(a))  a <- 0
#precision_mean_P <- c(3, 15, 100, 1000)[1 + a %% 4]
#mean_mean_P      <- c(3, 5, 7)[1 + a %/% 4]

# prior for MEAN_P is centered at 'mean_P_mean' and has strength 'mean_P_precision'
mean_P_precision <- 150
mean_P_mean      <- 7
lambda_P_mean <- (mean_P_mean / gamma(1+1/shape_P))^(-shape_P)


prior <- list(
  # flat priors
  #rate_P = 0.01, shape_P = 1,
  a_beta = 5, b_beta  = 5,  # beta(a_beta, b_beta) prior on beta
  
  # informative priors
  rate_P = (mean_P_precision - 1) / lambda_P_mean, shape_P = mean_P_precision, # gamma(shape_P, rate_P) prior on the Weibull (exponential) rate for P
  #a_beta = 1, b_beta = 1 # beta(a_beta, b_beta) prior on beta
  
  # other priors
  a_psi  = 1 , b_psi   = 1,   # beta(a_psi , b_psi ) prior on psi
  rate_H = 0.01, shape_H = 1 # gamma(shape_H, rate_H) prior on the Weibull rate for H
  
)

# initial values
mean_H_0 <- 130
rate_H_0 <- (mean_H_0/gamma(1+1/shape_H))^(-shape_H) # this gives a Weibull with mean mean_H_0
mean_P_0 <- mean_P_mean
rate_P_0 <- (mean_P_0/gamma(1+1/shape_P))^(-shape_P) # this gives a Weibull with mean mean_H_0
theta_0  <- list( # initial values for theta
  rate_H = rate_H_0, shape_H = shape_H,
  rate_P = rate_P_0, shape_P = shape_P,
  beta   = 0.8, psi = 0.1
) %>%
  update_scales()

#
# MCMC ####
out <- MCMC_cpp(
    d_obs_screen, d_obs_censor,
    theta_0, prior,
    epsilon_rate_H, epsilon_rate_P, epsilon_psi,
    t0, M, thin
    )
beep()


#
# Output ####
path_mcmc  <- paste0("output/MCMC/"   , data_origin)
path_fig   <- paste0("output/figures/", data_origin)
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
theta <- NULL
save(
    out,
    M, thin, n_cpu,
    theta, 
    prior, theta_0, epsilon_rate_H, epsilon_rate_P, epsilon_psi,
    d_obs_screen, d_obs_censor,
    file = file_draws
)

# save parameters only
theta_mcmc = out$THETA
save(
  theta_mcmc,
  file = paste0(path_mcmc, "/theta - ", sim_id)
)
