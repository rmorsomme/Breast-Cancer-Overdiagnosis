rm(list = ls())

path_mcmc   <- "output/MCMC"

setwd("C:/Users/18582/Desktop/Research/Marc/Breast Cancer Overdiagnosis/slurm")
n <- 42e3
m <- 1e2
factor_0 <- 1
N_CPU <- c(1, 2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
l <- length(N_CPU)
runtime <- numeric(l)

for(i in 1 : l){
    
    n_cpu <- N_CPU[i]
    
    file_id     <- paste0(path_mcmc, "/MCMC-n=", n, "-m=", m, "-factor0=", factor_0, "-n_cpu=", n_cpu, ".RDATA")
    load(file_id)
    
    runtime[i] <- out$runtime
    
}
plot(N_CPU, runtime, main = paste0("n=", n, ", m=", m))
