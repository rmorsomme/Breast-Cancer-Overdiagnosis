rm(list = ls())

path_mcmc   <- "output/MCMC"

setwd("C:/Users/18582/Desktop/Research/Marc/Breast Cancer Overdiagnosis/slurm")
n <- 10e3
M <- 5e3
factor_0 <- 1
n_cpu <- 5
EPSILON <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1)
l <- length(EPSILON)

rate_accept_psi <- numeric(l)
ESS_psi <- numeric(l)

for(i in 1 : l){ print(i)
    
    epsilon <- EPSILON[i]
    
    file_id <- paste0(
        path_mcmc, 
        "/MCMC-n=" , n, 
        "-M="      , M, 
        "-factor0=", factor_0, 
        "-n_cpu="  , n_cpu, 
        "-epsilon=", epsilon, 
        ".RDATA"
    )
    
    load(file_id)
    
    rate_accept_psi[i] <- mean(out$ACCEPT_PSI)
    ESS_psi[i]         <- coda::effectiveSize(out$THETA$PSI) 
    
}

plot(EPSILON, rate_accept_psi, main = paste0("n=", n, ", M=", M))
plot(EPSILON, ESS_psi        , main = paste0("n=", n, ", M=", M))
