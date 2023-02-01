
source("mcmc_setup.R")

out <- MCMC(d_process, d_obs_screen, d_obs_censor, theta_0, prior, m)

save(
    out,
    prior, theta,
    file = file_id
    )