
#
# setup ####
rm(list=ls())
set.seed(0)

library(tidyverse)
library(Rcpp)
library(coda)

source("slurm/R/helpers.r")

# helper

# computes the probability of dying in the interval [t0, t1].
compute_prob_death = function(t0, t1, OChaz) {
    
    t0-0.01 # add small value to 
    t1+0.01
    
    # endpoints of intervals
    if(t1 < ceiling(t0)) { # t1 is within the "same unit" as t0, e.g. t0=40.3 and t1=40.7
        ages = c(t0, t1)
    } else {
        ages = c(t0, ceiling(t0):floor(t1), t1)
    }
    
    # vector of probabilities of dying in each interval
    prob_death_interval = diff(ages) * OChaz$Estimate[ceiling(t0) : (floor(t1) + 1) + 1]
    
    # vector of probabilities of surviving in each interval
    prob_survival_interval = 1-prob_death_interval
    
    # probability of surviving all the intervals
    prob_survival = prod(prob_survival_interval)
    
    # probability of dying in any interval
    prob_death = 1 - prob_survival
    
    return(prob_death)
}

# load cause table
load("slurm/R/odx/OtherCause_birthcohort_1971.RDATA")

# paths to files and figures
path_mcmc  <- "slurm/output/MCMC/BCSC"
path_fig   <- "slurm/output/figures/BCSC"

# load MCMC draws
shape_P=2
sim_id     <- paste0(
    "M=", 1e3,
    "-AFS_low=", 50,
    "-AFS_upp=", 74,
    "-shape_H=", 2,
    "-shape_P=", shape_P,
    "-t0=", 30,
    ".RDATA"
)

file_draws <- paste(path_mcmc, sim_id, sep = "/")

load(file_draws)


# Take the subset of screen-detected individuals
ind_screen = which(d_obs_censor$censor_type=="screen")
n_screen   = length(ind_screen)
d_obs_censor_screen_age = d_obs_censor[ind_screen,]$censor_time

scale_P      = rate2scale(out$THETA$RATE_P, shape=shape_P)
Z_tau_screen = out$age_at_tau_hp_hat[,ind_screen]
Z_I_screen   = out$INDOLENT[,ind_screen]

# containers
S=length(scale_P)
prob_odx = tau_PC = matrix(nrow=S, ncol=n_screen)

for(i in 1:n_screen) { # loop over (screen-detected) individuals
    for(s in 1:M) { # loop over MCMC samples
        if(Z_I_screen[s,i]) {
            tau_PC[s,i] = Inf
            prob_odx[s,i] = 1
        } else {
            
            # TODO: generate multiple tau_PC for each (i,m)
            
            tau_PC[s,i] = Z_tau_screen[s,i] + rweibull(1, shape_P, scale_P[s]) # generate tau_PC such that tau_PC > censoring age
            while(tau_PC[s,i] < d_obs_censor_screen_age[i]) {
                tau_PC[s,i] = Z_tau_screen[s,i] + rweibull(1, shape_P, scale_P[s])
            }
            prob_odx[s,i] = compute_prob_death(t0=d_obs_censor_screen_age[i], t1 = tau_PC[s,i], OChaz)
        }
    }
}

# overall ODX rate
mean(prob_odx)

# have fun, be creative
prob_odx_ind = colMeans(prob_odx)

integers_age = which(round(d_obs_censor_screen_age) == d_obs_censor_screen_age)

mean(prob_odx_ind[-integers_age])
mean(prob_odx_ind[integers_age])

tibble(
    age_screen_detected = d_obs_censor_screen_age,
    prob_odx_ind = prob_odx_ind,
    prob_indolent = colMeans(Z_I_screen)
) %>%
    ggplot(aes(age_screen_detected, prob_odx_ind, col=prob_indolent)) +
    geom_point() +
    labs(title="ODX estimate per age of screen detection")
