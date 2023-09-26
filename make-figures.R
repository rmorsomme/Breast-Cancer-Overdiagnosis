
source("slurm/R/helpers.R")
source("slurm/R/helpers-figures.R")

library(tidyverse)
theme_set(theme_bw())

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


