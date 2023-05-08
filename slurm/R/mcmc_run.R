
{
    source("mcmc_setup.R")
    
    print(paste0("Number of CPUs: ", n_cpu))
    
    out <- MCMC(d_obs_screen, d_obs_censor, theta_0, prior, epsilon, M, thin, n_cpu, verbose = TRUE)
    
    beep()
}

save(
    out,
    prior, theta, theta_0,
    d_obs_screen, d_obs_censor,
    file = file_draws
    )
