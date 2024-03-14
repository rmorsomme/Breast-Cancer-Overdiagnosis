
#
# setup ####
rm(list=ls())

library(tidyverse)
library(mvtnorm)
library(Rcpp)
library(coda)

source("slurm/R/helpers.r")
source("slurm/R/helpers_for_calling_Rcpp.R")
sourceCpp("slurm/R/helpers_in_Rcpp.cpp")

results = tibble(
    shape_H=numeric(), shape_P=numeric(),
    lpd=numeric(), elpd_loo=numeric(), elpd_waic=numeric(), p_waic=numeric(),
    elpd_loo_vec=list()
)


#
## model setup ####
t0=30
for(shape_H in c(1,2)){
    for(shape_P in c(1,2)){  print(paste0("shape_H=",shape_H," - shape_P=",shape_P))
        
        #
        ## load MCMC output ####
        {
            # paths to files and figures
            path_mcmc  <- "slurm/output/MCMC/BCSC"
            path_fig   <- "slurm/output/figures/BCSC"
            
            # Load data
            sim_id     <- paste0(
                "M=", 1e3,
                "-AFS_low=", 50,
                "-AFS_upp=", 74,
                "-shape_H=", shape_H,
                "-shape_P=", shape_P,
                "-t0=", t0,
                ".RDATA"
            )
            
            file_draws <- paste(path_mcmc, sim_id, sep = "/")
            
            load(file_draws)
            
            # MCMC draws
            theta_MCMC = as_tibble(out$THETA)
            
            # make data object
            data.obj = make_dat.obj(d_obs_screen, d_obs_censor)
            
            # group indices
            screens  <- d_obs_censor$censor_type == "screen"
            censored <- d_obs_censor$censor_type == "censored"
            clinical <- d_obs_censor$censor_type == "clinical"
        }
        
        
        #
        # loo ####
        
        ## setup ####
        S = length(out$THETA$RATE_H)
        n = ncol(out$age_at_tau_hp_hat)
        
        ## compute loglik ####
        loglik = matrix(NA, nrow = S, ncol = n)
        for(s in 1:S){  #print(paste0(s, "/",S))
            # get theta
            theta_tbl = theta_MCMC[s,]
            theta_s = list(
                    rate_H = theta_tbl$RATE_H, shape_H = shape_H,
                    rate_P = theta_tbl$RATE_P, shape_P = shape_P,
                    beta   = theta_tbl$BETA, psi = theta_tbl$PSI
                ) %>% 
                update_scales()
            
            # compute loglik for each group
            # group: screen-detected
            age_at_tau_hp_hats_s = out$age_at_tau_hp_hat[s,screens]
            indolent_s = out$INDOLENT[s,screens]
            loglik[s,screens] = 
                dloglik_screens_obj(data.obj[[1]], theta_s, age_at_tau_hp_hats_s) +
                dloglik_sojourn_P_obj(data.obj[[1]], theta_s, age_at_tau_hp_hats_s, indolent_s)
            # group: censored
            age_at_tau_hp_hats_s = out$age_at_tau_hp_hat[s,censored]
            indolent_s = out$INDOLENT[s,censored]
            loglik[s,censored] = 
                dloglik_screens_obj(data.obj[[2]], theta_s, age_at_tau_hp_hats_s) +
                dloglik_sojourn_P_obj(data.obj[[2]], theta_s, age_at_tau_hp_hats_s, indolent_s)
            # group: interval-detected
            age_at_tau_hp_hats_s = out$age_at_tau_hp_hat[s,clinical]
            indolent_s = out$INDOLENT[s,clinical]
            loglik[s,clinical] = 
                dloglik_screens_obj(data.obj[[3]], theta_s, age_at_tau_hp_hats_s) +
                dloglik_sojourn_P_obj(data.obj[[3]], theta_s, age_at_tau_hp_hats_s, indolent_s)
        }
        
        ## compute accuracy measures ####
        
        ### lpd ####
        lik = exp(loglik)
        lpd = sum(log(colSums(lik)/S))
        
        ### elpd_loo ####
        lik_inv      = 1/lik
        p_loo        = S/colSums(lik_inv)
        elpd_loo_vec = log(p_loo)
        elpd_loo     = sum(elpd_loo_vec)
        
        ### WAIC ####
        V = numeric(n)
        for(i in 1:n)  V[i] = var(loglik[,i])
        p_waic = sum(V)
        elpd_waic = lpd - p_waic
        
        ## summarize results ####
        results = results %>%
            add_row(
                shape_H=shape_H, shape_P=shape_P,
                lpd=lpd, elpd_loo=elpd_loo, elpd_waic=elpd_waic, p_waic=p_waic, 
                elpd_loo_vec=list(elpd_loo_vec)
            )
    }
}

# compare the elpd_loo with model 4
se = numeric(3)
for(model in 1:3) {
    se    [model] = sqrt(n*var(results$elpd_loo_vec[[4]] - results$elpd_loo_vec[[model]]))
}

# t statistics
t_stat = abs(results$elpd_loo[4] - results$elpd_loo[1:3]) / se
# p-values
p_val  = pt(t_stat, df = n-1, lower.tail = FALSE)

results = results %>%
    mutate(
        se = c(se, NA),
        t_stat = c(t_stat, NA),
        p_val = c(p_val, NA),
    )
results



