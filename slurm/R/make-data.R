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
shape_H <- 2
shape_P <- 1


#
# Parameters ####
# these values are taken from Ryser et al. (2022)
theta <- list(
    rate_H = 7.5e-5, shape_H = shape_H,
    rate_P = 1/7  , shape_P = shape_P,
    beta   = 38.5/(38.5+5.8), psi = 0.1
) %>% update_scales()


#
## Generate data ####
n  <- 37e3 #37e3  # sample size
source("R/simulator.R")


#
# Save data ####
path_data  <- "output/MCMC/simulation/data"
id_data     <- paste0(
    "n=", n,
    "-shape_H=", shape_H,
    "-shape_P=", shape_P,
    "-t0=", t0,
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
