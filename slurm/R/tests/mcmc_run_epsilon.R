
source("mcmc_setup.R")

EPSILON <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1)

epsilon_i <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
epsilon   <- EPSILON[epsilon_i]

print(epsilon)

out <- MCMC(d_process, d_obs_screen, d_obs_censor, theta_0, prior, epsilon, M, thin, n_cpu)

file_id <- paste0(
    path_mcmc, 
    "/MCMC-n=" , n, 
    "-M="      , M, 
    "-factor0=", factor_0, 
    "-n_cpu="  , n_cpu, 
    "-epsilon=", epsilon, 
    ".RDATA"
    )

save(
    out,
    prior, theta,
    d_process, d_obs_screen, d_obs_censor,
    file = file_id
)
