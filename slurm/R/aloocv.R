
#
# setup ####
rm(list=ls())

set.seed(0)

library(tidyverse)
library(mvtnorm)
library(Rcpp)
library(coda)


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
    data.obj = make_dat.obj(d_obs_screen %>% filter(person_id==105), d_obs_censor %>% filter(person_id==105))
    
    # group indices
    screens  <- d_obs_censor$censor_type == "screen"
    censored <- d_obs_censor$censor_type == "censored"
    clinical <- d_obs_censor$censor_type == "clinical"
    
    # unique person id
    id_vec = unique(d_obs_censor$person_id) %>% sample()
}


#
# loo ####

## setup ####
S = nrow(theta_MCMC)
n = nrow(d_obs_censor)
n=1e4
J=200 # number of IS samples
ess_target = 100
# current: n=10e3 and J=2e3 takes 17 hours per model
# goal: n=40e3 and J=3e3
z_tau_HP=z_I=qlog=loglik_joint=numeric(J) # containers for IS

# lik 
lik = ess = n_is = matrix(0, nrow = S, ncol = n)

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
    
    # get data for subject i
    id = id_vec[i]
    group = d_obs_censor %>% filter(person_id==id) %>% pull(censor_type)
    
    d_obs_censor_i = d_obs_censor %>% filter(person_id==id)
    d_obs_screen_i = d_obs_screen %>% filter(person_id==id)
    
    # replicate the data J times
    d_obs_censor_i_rep = tibble(
        person_id = 1:J,
        AFS = d_obs_censor_i$AFS,
        censor_type = d_obs_censor_i$censor_type,
        censor_time = d_obs_censor_i$censor_time
    )
    
    d_obs_screen_i_rep = tibble(
        person_id = rep(1:J, each=nrow(d_obs_screen_i)),
        AFS = rep(d_obs_screen_i$AFS, J),
        screen_id = rep(d_obs_screen_i$screen_id, J),
        age_screen = rep(d_obs_screen_i$age_screen, J),
        screen_detected = rep(d_obs_screen_i$screen_detected, J)
    )
    
    data_object = make_dat.obj(d_obs_screen_i_rep, d_obs_censor_i_rep) %>% .[[group]]
    
    for(s in 1:S) {  #print(paste0(s, "/",S))
        
        integrand = NULL
        n_IS = 0
        
        while(ess[s,i] < ess_target && n_IS < 25) {
            
            # importance sampling
            prob_tau      = compute_prob_tau_obj(data_object, theta_list[[s]], t0)
            z_tau_HP      = rprop_age_at_tau_hp_hat_obj(data_object, prob_tau, theta_list[[s]], t0)
            prob_indolent = compute_prob_indolent_obj(data_object, theta_list[[s]], z_tau_HP)
            z_I           = rprop_indolent_obj(data_object, prob_indolent)
            qlog          = dlog_prop_latent_obj(data_object, prob_tau, theta_list[[s]], z_tau_HP, prob_indolent, z_I, t0)
            loglik_joint  = dlog_likelihood_obj(data_object, theta_list[[s]], z_tau_HP, z_I, t0)
            
            # IS estimate
            integrand = c(integrand, exp(loglik_joint-qlog))
            lik[s,i]  = mean(integrand)
            
            # number of IS samples needed
            n_IS = n_IS+1
            
            # ESS
            integrand_norm=integrand / sum(integrand)
            ess[s,i]=1/sum(integrand_norm^2)
            
        }
        n_is[s,i] = n_IS * J
    }
}


results = tibble(
    shape_H=shape_H, shape_P=shape_P,
    lik=list(lik), ess=list(ess), n_is=list(n_is)
)

save(results, file=paste0(path_aloocv, "/shape_H=",shape_H,"-shape_P=",shape_P,".RDATA"))

        # for(s in 1:S){  
        #     
        #     # get theta
        #     theta_tbl = theta_MCMC[s,]
        #     theta_s = list(
        #         rate_H = theta_tbl$RATE_H, shape_H = shape_H,
        #         rate_P = theta_tbl$RATE_P, shape_P = shape_P,
        #         beta   = theta_tbl$BETA, psi = theta_tbl$PSI
        #     ) %>% 
        #         update_scales()
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
