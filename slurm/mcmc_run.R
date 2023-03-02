
{
    source("mcmc_setup.R")
    
    print(n_cpu)
    
    out <- MCMC(d_process, d_obs_screen, d_obs_censor, theta_0, prior, epsilon, m, n_cpu = n_cpu)
    
    beep()
}

save(
    out,
    prior, theta,
    d_process, d_obs_screen, d_obs_censor,
    file = file_id
    )