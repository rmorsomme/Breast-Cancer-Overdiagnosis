
source("mcmc_setup.R")

N_CPU <- c(1, 2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)

for(n_cpu in N_CPU){  print(n_cpu)
    
    out <- MCMC(d_process, d_obs_screen, d_obs_censor, theta_0, prior, epsilon, M, thin, n_cpu)
    
    file_id <- paste0(path_mcmc, "/MCMC-n=", n, "-M=", M, "-factor0=", factor_0, "-n_cpu=", n_cpu, ".RDATA")
    
    save(
        out,
        prior, theta,
        d_process, d_obs_screen, d_obs_censor,
        file = file_id
    )
    
}
