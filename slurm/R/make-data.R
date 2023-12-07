#
# Set up
#setwd("C:/Users/18582/Desktop/Research/Marc/Breast Cancer Overdiagnosis/slurm")
rm(list = ls())

library(tidyverse)
library(tictoc)
library(beepr)
library(Rcpp)

sourceCpp("R/helpers_in_Rcpp.cpp")
source("R/helpers_for_calling_Rcpp.R")

set.seed(1)

#
# Model setup ####
t0 <- 30
parameter_from_Ryser2022 = FALSE

#
# Parameters ####
if(parameter_from_Ryser2022) {
    
    # these values are taken from Ryser et al. (2022)
    theta <- list(
        rate_H = 7.5e-5, shape_H = 2,
        rate_P = 1/7  , shape_P = 1,
        beta   = 38.5/(38.5+5.8), psi = 0.1
    ) %>% update_scales()
    
} else {
    
    shape_H <- 2
    shape_P <- 2
    
    beta = 0.8
    psi  = 0.5
    
    mean_H <- 100 # expected number of years
    rate_H <- (mean_H/gamma(1+1/shape_H))^(-shape_H) # this gives a Weibull with mean mean_H_0
    mean_P <- 10 # expected number of years
    rate_P <- (mean_P/gamma(1+1/shape_P))^(-shape_P) # this gives a Weibull with mean mean_H_0
    
    theta <- list(
        rate_H = rate_H, shape_H = shape_H,
        rate_P = rate_P, shape_P = shape_P,
        beta   = beta  , psi = psi
    ) %>% update_scales()
}


#
## Generate data ####
n  <- 5e4 #37e3  # sample size
source("R/simulator.R")


#
# Save data ####
path_data <- "data/simulated"
id_data   <- paste0(
    "n=", n,
    "-shape_H=", shape_H,
    "-shape_P=", shape_P,
    "-t0=", t0,
    "-psi=", psi,
    "-mean_P=", mean_P,
    ".RDATA"
)

file_data <- paste(path_data, id_data, sep = "/")
save(
    n, shape_H, shape_P, t0,
    theta,
    d_obs_screen, d_obs_censor,
    file = file_data
)

beep()
