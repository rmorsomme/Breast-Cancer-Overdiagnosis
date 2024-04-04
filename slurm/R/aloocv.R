
#
# setup ####
rm(list=ls())

library(tidyverse)
library(Rcpp)


a <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if(is.na(a)){
    a = 3
    setwd("C:/Users/18582/Desktop/Research/Marc/Breast Cancer Overdiagnosis/slurm")
}
shape_H = c(1,2)[1+a%%2]
shape_P = c(1,2)[1+a%/%2]
print(paste0("shape_H=",shape_H," - shape_P=",shape_P))

source("R/helpers.R")
source("R/helpers_for_calling_Rcpp.R")
sourceCpp("R/helpers_in_Rcpp.cpp")

path_aloocv <- "output/misc/aloocv"


#
## load MCMC output ####
{
    # paths to files and figures
    path_mcmc  <- "output/MCMC/BCSC"
    path_fig   <- "output/figures/BCSC"
    
    # Load data
    t0=30
    sim_id     <- paste0(
        "M=", 1e5,
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
    
    # unique person id
    id_vec = unique(d_obs_censor$person_id)
}


#
# loo ####

## setup ####
S = nrow(theta_MCMC)
n = nrow(d_obs_censor)
J_increment = 225 # number of IS samples per increment
ess_target = 200 # target ESS for each (i,s)


# containers
z_tau_HP = z_I = qlog = loglik_joint = numeric(J_increment)
lik = ess = J = matrix(0, nrow = S, ncol = n)

# prepare theta
theta_list = vector("list", S)
for(s in 1:S) {  #print(paste0(s, "/",S))
    
    # get theta
    theta_tbl = theta_MCMC[s,]
    theta_list[[s]] = list(
        rate_H = theta_tbl$RATE_H, shape_H = shape_H,
        rate_P = theta_tbl$RATE_P, shape_P = shape_P,
        beta   = theta_tbl$BETA, psi = theta_tbl$PSI
    ) %>% 
        update_scales()
}


# mcmapply(
#     dloglik_sojourn_H_i,
#     tau_HP, censor_time,
#     MoreArgs = list(theta = theta, t0 = t0),
#     USE.NAMES = FALSE,
#     mc.cores = n_cpu
# ) %>%
#     sum()
#SBATCH --nodes=n_cpu



## aloocv ####
for(i in 1:n) {  print(paste0(i, "/",n)) # to be parallelized
    
    # id of i-th individual
    id = id_vec[i]
    
    # get data for subject i
    group = d_obs_censor %>% filter(person_id==id) %>% pull(censor_type)
    
    d_obs_censor_i = d_obs_censor %>% filter(person_id==id)
    d_obs_screen_i = d_obs_screen %>% filter(person_id==id)
    
    # replicate the data J_increment times
    d_obs_censor_i_rep = tibble(
        person_id = 1:J_increment,
        AFS = d_obs_censor_i$AFS,
        censor_type = d_obs_censor_i$censor_type,
        censor_time = d_obs_censor_i$censor_time
    )
    
    d_obs_screen_i_rep = tibble(
        person_id = rep(1:J_increment, each=nrow(d_obs_screen_i)),
        AFS = rep(d_obs_screen_i$AFS, J_increment),
        screen_id = rep(d_obs_screen_i$screen_id, J_increment),
        age_screen = rep(d_obs_screen_i$age_screen, J_increment),
        screen_detected = rep(d_obs_screen_i$screen_detected, J_increment)
    )
    
    data_object = make_dat.obj(d_obs_screen_i_rep, d_obs_censor_i_rep) %>% .[[group]]
    
    for(s in 1:S) {  #print(paste0(s, "/",S))
        
        set.seed(s)
        
        integrand = NULL
        n_increment = 0
        
        while(ess[s,i] < ess_target && n_increment < 25) {
            
            # importance sampling
            prob_tau      = compute_prob_tau_obj(data_object, theta_list[[s]], t0)
            z_tau_HP      = rprop_age_at_tau_hp_hat_obj(data_object, prob_tau, theta_list[[s]], t0)
            prob_indolent = compute_prob_indolent_obj(data_object, theta_list[[s]], z_tau_HP)
            z_I           = rprop_indolent_obj(data_object, prob_indolent)
            qlog          = dlog_prop_latent_obj(data_object, prob_tau, theta_list[[s]], z_tau_HP, prob_indolent, z_I, t0)
            loglik_joint  = dlog_likelihood_obj(data_object, theta_list[[s]], z_tau_HP, z_I, t0)
            
            # IS estimate
            integrand = c(integrand, exp(loglik_joint - qlog))
            lik[s,i]  = mean(integrand)
            
            # ESS
            integrand_norm=integrand / sum(integrand)
            ess[s,i]=1/sum(integrand_norm^2)
            
            # number of increments
            n_increment = n_increment+1
            
        }
        J[s,i] = n_increment * J_increment
    }
}


results = tibble(
    shape_H=shape_H, shape_P=shape_P, lik=list(lik), ess=list(ess), J=list(J)
    )
results_summary = results %>%
    mutate(
        ess_min = ess %>% map_dbl(min),
        ess_mean = ess %>% map_dbl(mean),
        ess_q001 = ess %>% map_dbl(quantile, probs=0.01),
        J_max = J %>% map_dbl(min),
        J_mean = J %>% map_dbl(mean),
        J_q099 = J %>% map_dbl(quantile, probs=0.99),
        lpd = lik %>% map_dbl(~sum(log(colMeans(.)))),
        elpd_loo_vec = lik %>% map(~log(1/colMeans(1/.))),
        elpd_loo = elpd_loo_vec %>% map_dbl(mean)
    ) %>%
    select(-lik, -ess, -J)

save(results        , file=paste0(path_aloocv, "/shape_H=",shape_H,"-shape_P=",shape_P,".RDATA"))
save(results_summary, file=paste0(path_aloocv, "/summary-shape_H=",shape_H,"-shape_P=",shape_P,".RDATA"))

        # for(s in 1:S){  
        #     
        #     
        #     # compute loglik for each group
        #     # group: screen-detected
        #     age_at_tau_hp_hats_s = out$age_at_tau_hp_hat[s,screens]
        #     indolent_s = out$INDOLENT[s,screens]
        #     loglik[s,screens] = 
        #         dloglik_screens_obj(data.obj[[1]], theta_s, age_at_tau_hp_hats_s) +
        #         dloglik_sojourn_P_obj(data.obj[[1]], theta_s, age_at_tau_hp_hats_s, indolent_s)
        #     # group: censored
        #     age_at_tau_hp_hats_s = out$age_at_tau_hp_hat[s,censored]
        #     indolent_s = out$INDOLENT[s,censored]
        #     loglik[s,censored] = 
        #         dloglik_screens_obj(data.obj[[2]], theta_s, age_at_tau_hp_hats_s) +
        #         dloglik_sojourn_P_obj(data.obj[[2]], theta_s, age_at_tau_hp_hats_s, indolent_s)
        #     # group: interval-detected
        #     age_at_tau_hp_hats_s = out$age_at_tau_hp_hat[s,clinical]
        #     indolent_s = out$INDOLENT[s,clinical]
        #     loglik[s,clinical] = 
        #         dloglik_screens_obj(data.obj[[3]], theta_s, age_at_tau_hp_hats_s) +
        #         dloglik_sojourn_P_obj(data.obj[[3]], theta_s, age_at_tau_hp_hats_s, indolent_s)
        # }
