setwd("C:/Users/18582/Desktop/Research/Marc/Breast Cancer Overdiagnosis/slurm")
library(tidyverse)
rm(list = ls())

M     <- 2e3
n_cpu <- 15
path_mcmc  <- "output/MCMC/BCSC"
path_fig   <- "output/figures/BCSC"
shape_H <- 2 # 1, 1.1, 1.5, 2 ,2.5, 3, 3.5
shape_P <- 1

for(precision_mean_P in c(3, 15, 100, 1000)){ print(precision_mean_P)
    for(mean_mean_p in c(3, 5, 7)){
        for(AFS in 45){
            for(t0 in c(30, 40, AFS - 1)){
                sim_id     <- paste0(
                    "M=", M,
                    "-AFS=", AFS,
                    "-shape_H=", shape_H, # 1, 1.1, 1.5, 2 ,2.5, 3, 3.5
                    "-shape_P=", shape_P,
                    "-t0=", t0,
                    "-mean_mean_p=", mean_mean_p,
                    "-precision_mean_P=", precision_mean_P,
                    "-n_cpu=", n_cpu, 
                    ".RDATA"
                )
                file_draws <- paste(path_mcmc, sim_id, sep = "/")
                file_fig   <- paste(path_fig , sim_id, sep = "/")
                setwd("C:/Users/18582/Desktop/Research/Marc/Breast Cancer Overdiagnosis/slurm")
                load(file_draws)
                
                theta$mean_H <- theta$rate_H^(-1/shape_H)*gamma(1+1/shape_H)
                theta$mean_P <- theta$rate_P^(-1/shape_P)*gamma(1+1/shape_P)
                
                THETA <- out$THETA %>%
                    as_tibble() %>%
                    mutate(
                        MEAN_H = RATE_H^(-1/theta_0$shape_H)*gamma(1+1/theta_0$shape_H),
                        MEAN_P = RATE_P^(-1/theta_0$shape_P)*gamma(1+1/theta_0$shape_P),
                        iteration = 1:nrow(.)
                    ) %>%
                    filter(iteration > 100) %>% # remove burnin
                    pivot_longer(cols = -iteration, names_to = "parameter", values_to = "draws")
                THETA_summary <- THETA %>%
                    group_by(parameter) %>%
                    summarize(
                        q_low   = quantile(draws, probs = 0.05),
                        q_high  = quantile(draws, probs = 0.95)
                    )
                
                ggplot(THETA, aes(draws)) +
                    geom_histogram(aes(y = after_stat(density))) +
                    facet_wrap(~parameter, scales = "free") +
                    geom_vline(data = THETA_summary, mapping = aes(xintercept=q_low)) +
                    geom_vline(data = THETA_summary, mapping = aes(xintercept=q_high)) +
                    labs(
                        title = paste0("AFS=", AFS, " - t0=", t0),
                        subtitle = paste0("Prior on `mean_p`: mean=", mean_mean_p,", strength=", precision_mean_P)
                    )
                ggsave(
                    paste(sim_id, "histogram.jpg", sep = "_"),
                    path = path_fig, width = 1.61803, height = 1, scale = 5
                )
            }
        }
    }
}


