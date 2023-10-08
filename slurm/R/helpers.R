
#
# Weibull ####

rate2scale <- function(rate, shape){
    rate^(-1/shape)
}

update_scales <- function(theta){
    theta$scale_H <- rate2scale(theta$rate_H, theta$shape_H)
    theta$scale_P <- rate2scale(theta$rate_P, theta$shape_P)
    return(theta)
}

pweibull_ab <- function(a, b, shape, scale){
    # exp(-(a/scale)^shape) -  exp(-(b/scale)^shape) (slower than pweibull())
    pweibull(b, shape, scale, lower.tail = T) -
        pweibull(a, shape, scale, lower.tail = T)
}

rweibull_trunc <- function(a, b, shape, scale){
    u <- runif(1, pweibull(a, shape, scale), pweibull(b, shape, scale))
    qweibull(u, shape, scale)
}

dweibull_trunc <- function(x, a, b, shape, scale, log = T){
    if(log){
        dweibull(x, shape, scale, log = T) - log(pweibull_ab(a, b, shape, scale))
    }else{
        dweibull(x, shape, scale, log = F) /     pweibull_ab(a, b, shape, scale)
    }
}

#
# theta ####

add_rate_H <- function(theta, rate_H){
    theta$rate_H <- rate_H
    theta <- update_scales(theta)
    return(theta)
}

add_rate_P <- function(theta, rate_P){
    theta$rate_P <- rate_P
    theta <- update_scales(theta)
    return(theta)
}

add_beta <- function(theta, beta){
    theta$beta <- beta
    return(theta)
}

add_psi <- function(theta, psi){
    theta$psi <- psi
    return(theta)
}

#
# Gibbs theta ####

gibbs_rate_H <- function(tau_HP, theta, prior, censor_time, t0){ # gamma prior and weibull likelihood conjugacy
    
    a_n <- prior$shape_H + sum(is.finite(tau_HP)) # number of observed weibull
    
    sojourn_H   <- pmin(tau_HP, censor_time) - t0
    sum_sojourn <- sum(sojourn_H^theta$shape_H) 
    b_n         <- prior$rate_H + sum_sojourn
    
    rgamma(1, shape = a_n, rate = b_n)
}

gibbs_rate_P <- function(tau_HP, indolent, theta, prior, censor_type, censor_time){ # gamma prior and weibull likelihood conjugacy
    
    a_n <- prior$shape_P + sum(censor_type == "clinical") # number of observed weibull
    
    sojourn_P   <- (censor_time - tau_HP)[!indolent & is.finite(tau_HP)]
    sum_sojourn <- sum(sojourn_P^theta$shape_P) 
    b_n         <- prior$rate_P + sum_sojourn
    
    rgamma(1, shape = a_n, rate = b_n)
}

gibbs_psi <- function(indolent, prior){ # beta prior and binomial likelihood conjugacy
    
    a_n <- prior$a_psi + sum(indolent  , na.rm = TRUE) 
    b_n <- prior$b_psi + sum(1-indolent, na.rm = TRUE) 
    
    rbeta(1, a_n, b_n)
}

gibbs_beta <- function(tau_HP, theta, age_screen, prior, n_screen_positive){ # beta prior and binomial likelihood conjugacy
    
    n_screen <- list(tau_HP, age_screen) %>% # TODO: parallelize
        pmap_dbl(~sum(..1 < ..2)) %>%
        sum()
    
    a_n <- prior$a_beta + n_screen_positive
    b_n <- prior$b_beta + (n_screen-n_screen_positive)
    
    beta_new <- rbeta(1, a_n, b_n)
    
    theta <- add_beta(theta, beta_new)
    return(theta)
    
}

#
# Likelihood ####

dloglik_sojourn_H_i <- function(tau_HP_i, theta, censor_time_i, t0){
    
    if(is.infinite(tau_HP_i)){
        pweibull(
            censor_time_i - t0,
            shape = theta$shape_H, scale = theta$scale_H,
            log.p = TRUE, lower.tail = FALSE
        )
    }else if(is.finite(tau_HP_i)){
        dweibull(
            tau_HP_i - t0, 
            shape = theta$shape_H, scale = theta$scale_H,
            log = TRUE
        )
    }
    
}

dloglik_sojourn_H <- function(tau_HP, theta, censor_time, t0, n_cpu){
    
    mcmapply(
        dloglik_sojourn_H_i,
        tau_HP, censor_time,
        MoreArgs = list(theta = theta, t0 = t0),
        USE.NAMES = FALSE,
        mc.cores = n_cpu
    ) %>%
        sum()
    
}

dloglik_sojourn_P_i <- function(tau_HP_i, theta, indolent_i, censor_type_i, censor_time_i){
    
    if(is.infinite(tau_HP_i) | indolent_i){
        0
    }else if(censor_type_i == "clinical"){
        dweibull(
            censor_time_i - tau_HP_i, 
            shape = theta$shape_P, scale = theta$scale_P, 
            log = TRUE
        )
    } else if(censor_type_i %in% c("screen", "censored")){
        pweibull(
            censor_time_i - tau_HP_i, 
            shape = theta$shape_P, scale = theta$scale_P, 
            log.p = TRUE, lower.tail = FALSE
        )
    }
    
}

dloglik_sojourn_P <- function(tau_HP, theta, indolent, censor_type, censor_time, n_cpu){
    
    mcmapply(
        dloglik_sojourn_P_i,
        tau_HP, indolent, censor_type, censor_time,
        MoreArgs = list(theta = theta),
        USE.NAMES = FALSE,
        mc.cores = n_cpu
    ) %>%
        sum()
    
}

dloglik_indolent <- function(indolent_i, theta){
    if(indolent_i) { # Bernoulli pmf
        log(theta$psi)
    } else if(!indolent_i){
        log(1-theta$psi)
    }
}

dloglik_screens <- function(age_screen_i, theta, tau_HP_i, n_screen_positive_i){
    
    n_screen_i          <- sum(age_screen_i > tau_HP_i) # number of screens during pre-clinical phase
    n_screen_negative_i <- n_screen_i - n_screen_positive_i
    
    # product of Bernoulli RVs
    n_screen_positive_i*log(theta$beta) + n_screen_negative_i*log(1-theta$beta)
    
}

dlog_likelihood_i <- function(tau_HP_i, indolent_i, censor_type_i, censor_time_i, age_screen_i, n_screen_positive_i, theta, t0){
    
    dlog_H <- dloglik_sojourn_H_i(tau_HP_i, theta, censor_time_i, t0)
    dlog_P <- dloglik_sojourn_P_i(tau_HP_i, theta, indolent_i, censor_type_i, censor_time_i)
    dlog_S <- dloglik_screens (age_screen_i, theta, tau_HP_i, n_screen_positive_i)
    dlog_I <- dloglik_indolent(indolent_i, theta)
    
    dlog_lik <- dlog_H + dlog_P + dlog_S + dlog_I
    return(dlog_lik)
    
}


dlog_likelihood <- function(tau_HP, indolent, censor_type, censor_time, age_screen, n_screen_positive, theta, n_cpu, t0){
    mcmapply(
        dlog_likelihood_i,
        tau_HP, indolent, censor_type, censor_time, age_screen, n_screen_positive,
        MoreArgs = list(theta = theta, t0 = t0),
        mc.cores = n_cpu
    ) %>% sum()
}

dloglik_PI_i <- function(tau_HP_i, indolent_i, censor_type_i, censor_time_i, theta){
    
    dlog_P <- dloglik_sojourn_P_i(tau_HP_i, theta, indolent_i, censor_type_i, censor_time_i)
    dlog_I <- dloglik_indolent   (indolent_i, theta)
    
    dlog_lik <- dlog_P + dlog_I
    return(dlog_lik)
    
}

dloglik_PI <- function(tau_HP, indolent, censor_type, censor_time, theta, n_cpu){
    mcmapply(
        dloglik_PI_i,
        tau_HP, indolent, censor_type, censor_time,
        MoreArgs = list(theta = theta),
        mc.cores = n_cpu
    ) %>% sum()
}


dloglik_psi <- function(tau_HP, indolent, censor_type, censor_time, theta, AFS, n_AFS, t0, n_cpu){
    
    dlog_cp <- dloglik_cp(AFS, n_AFS, theta, t0)
    dlog_PI <- dloglik_PI(tau_HP, indolent, censor_type, censor_time, theta, n_cpu)
    
    dlog_psi <- dlog_cp + dlog_PI
    return(dlog_psi)
    
}



dloglik_tau_i <- function(tau_HP_i, indolent_i, censor_type_i, censor_time_i, age_screen_i, n_screen_positive_i, theta, t0){
    
    dlog_H <- dloglik_sojourn_H_i(tau_HP_i, theta, censor_time_i, t0)
    dlog_P <- dloglik_sojourn_P_i(tau_HP_i, theta, indolent_i, censor_type_i, censor_time_i)
    dlog_S <- dloglik_screens(age_screen_i, theta, tau_HP_i, n_screen_positive_i)
    
    dlog_lik <- dlog_H + dlog_P + dlog_S
    return(dlog_lik)
    
}


dloglik_tau <- function(tau_HP, indolent, censor_type, censor_time, age_screen, n_screen_positive, theta, n_cpu, t0){
    mcmapply(
        dloglik_tau_i,
        tau_HP, indolent, censor_type, censor_time, age_screen, n_screen_positive,
        MoreArgs = list(theta = theta, t0 = t0),
        mc.cores = n_cpu
    ) %>% sum()
}


dloglik_cp <- function(AFS, n_AFS, theta, t0){
    
    cp_log_AFS <- AFS %>% map_dbl(compute_cp_log, theta = theta, t0 = t0) # TODO: parallelize
    
    dlog_cp <- n_AFS %*% (-cp_log_AFS)
    return(as.numeric(dlog_cp))
    
}



dloglik_rate_H <- function(tau_HP, theta, censor_time, AFS, n_AFS, t0, n_cpu){
    
    dlog_cp <- dloglik_cp       (AFS, n_AFS, theta, t0)
    dlog_H  <- dloglik_sojourn_H(tau_HP, theta, censor_time, t0, n_cpu)
    
    dlog_rate_H <- dlog_cp + dlog_H
    
    return(dlog_rate_H)
    
}


dloglik_rate_P <- function(tau_HP, theta, indolent, censor_type, censor_time, AFS, n_AFS, t0, n_cpu){
    
    dlog_cp <- dloglik_cp       (AFS, n_AFS, theta, t0)
    dlog_P  <- dloglik_sojourn_P(tau_HP, theta, indolent, censor_type, censor_time, n_cpu)
    
    dlog_rate_P <- dlog_cp + dlog_P
    
    return(dlog_rate_P)
    
}


#
# Conditioning probability ####

compute_integral <- function(L, U, theta){
    
    integrand <- function(t, U, shape_H, shape_P, scale_H, scale_P){
        dweibull(    t, shape_H, scale_H                    ) * 
        pweibull(U - t, shape_P, scale_P, lower.tail = FALSE)
    }
    
    integrate(
        integrand, L, U,
        U=U, 
        shape_H=theta$shape_H, shape_P=theta$shape_P, 
        scale_H=theta$scale_H, scale_P=theta$scale_P, 
        stop.on.error = TRUE
    ) %>% 
        .[["value"]]
}

compute_cp_log <- function(theta, AFS, t0){
    
    L <- 0        # lower bound
    U <- AFS - t0 # upper bound
    
    prob_onset_after <- pweibull(U, theta$shape_H, theta$scale_H, lower.tail = FALSE)
    
    integral <- compute_integral(L, U, theta)
    
    cp     <- theta$psi + (1-theta$psi) * (prob_onset_after + integral) # conditioning probability
    cp_log <- log(cp)
    return(cp_log)
    
}


#
# M-H rate_H ####

rprop_rate_H <- function(theta, epsilon_rate_H){
    
    rate_H <- theta$rate_H + runif(1, - epsilon_rate_H, epsilon_rate_H)
    rate_H <- abs(rate_H)
    
    return(rate_H)
    
}

MH_rate_H <- function(tau_HP, theta_cur, prior, censor_time, t0, AFS, n_AFS, n_cpu, epsilon_rate_H){
    
    # current value
    rate_H_cur  <- theta_cur$rate_H # for evaluating prior density
    
    # propose a new value
    rate_H_new <- rprop_rate_H(theta_cur, epsilon_rate_H) # symmetric proposal
    theta_new  <- add_rate_H(theta_cur, rate_H_new)  # for dloglik_rate_H()
    
    # M-H acceptance ratio
    dlog_lik_cur   <- dloglik_rate_H(tau_HP, theta_cur, censor_time, AFS, n_AFS, t0, n_cpu)
    dlog_lik_new   <- dloglik_rate_H(tau_HP, theta_new, censor_time, AFS, n_AFS, t0, n_cpu)
    dlog_prior_new <- dgamma(rate_H_new, shape = prior$shape_H, rate = prior$rate_H, log = TRUE)
    dlog_prior_cur <- dgamma(rate_H_cur, shape = prior$shape_H, rate = prior$rate_H, log = TRUE)
    
    MH_logratio <- (dlog_lik_new + dlog_prior_new) - (dlog_lik_cur + dlog_prior_cur)
    
    out <- if(runif(1) < exp(MH_logratio)){ # accept new value
        list(theta = theta_new, accept = TRUE)
    }else{ # keep current value
        list(theta = theta_cur, accept = FALSE)
    }
    
    return(out)
    
}


#
# M-H rate_P ####

rprop_rate_P <- function(theta, epsilon_rate_P){
    
    rate_P <- runif(1, theta$rate_P - epsilon_rate_P, theta$rate_P + epsilon_rate_P)
    rate_P <- abs(rate_P)
    
    return(rate_P)
    
}

MH_rate_P <- function(tau_HP, theta_cur, prior, indolent, censor_type, censor_time, t0, AFS, n_AFS, n_cpu, epsilon_rate_H){
    
    # current value
    rate_P_cur  <- theta_cur$rate_P # for dgamma()
    
    # propose a new value
    rate_P_new <- rprop_rate_P(theta_cur, epsilon_rate_P) # symmetric proposal on positive line
    theta_new  <- add_rate_P(theta_cur, rate_P_new)  # for dloglik_rate_P()
    
    # M-H acceptance ratio
    dlog_lik_cur   <- dloglik_rate_P(tau_HP, theta_cur, indolent, censor_type, censor_time, AFS, n_AFS, t0, n_cpu)
    dlog_lik_new   <- dloglik_rate_P(tau_HP, theta_new, indolent, censor_type, censor_time, AFS, n_AFS, t0, n_cpu)
    dlog_prior_new <- dgamma(rate_P_new, shape = prior$shape_P, rate = prior$rate_P, log = TRUE)
    dlog_prior_cur <- dgamma(rate_P_cur, shape = prior$shape_P, rate = prior$rate_P, log = TRUE)
    
    MH_logratio <- (dlog_lik_new + dlog_prior_new) - (dlog_lik_cur + dlog_prior_cur)
    
    out <- if(runif(1) < exp(MH_logratio)){ # accept new value
        list(theta = theta_new, accept = TRUE)
    }else{ # keep current value
        list(theta = theta_cur, accept = FALSE)
    }
    
    return(out)
    
}


#
# M-H tau_HP ####

compute_endpoints_i <- function(age_screen, censor_type_i, censor_time_i, t0){
    
    # construct endpoints of intervals (t_i, t_{i+1}]
    endpoints <- if(censor_type_i == "screen"){ # if screen, age_screen are endpoints
        c(t0, age_screen)
    }else if(censor_type_i == "clinical"){ # if clinical, add tau_PC as endpoint of last interval
        c(t0, age_screen, censor_time_i     )
    }else if(censor_type_i == "censored"){ # if censored, add tau_PC and Inf as endpoints
        c(t0, age_screen, censor_time_i, Inf)
    }
    
    return(endpoints)
    
}

compute_endpoints <- function(age_screen, censor_type, censor_time, n_cpu, t0){
    mcmapply(
        compute_endpoints_i,
        age_screen, censor_type, censor_time,
        MoreArgs = list(t0 = t0),
        mc.cores = n_cpu
    )
}

compute_prob_tau_i <- function(censor_type_i, censor_time_i, endpoints_i, theta, t0){
    
    K <- length(endpoints_i) - 1 # number of intervals
    
    # compute probability of interval based on the weibull waiting time in H (screen, censored) or in P (clinical)
    prob_interval <- if(censor_type_i == "screen"){
        
        pweibull_ab(
            endpoints_i[1:K    ] - t0,
            endpoints_i[1:K + 1] - t0,
            shape = theta$shape_H, scale = theta$scale_H
        )
        
    } else if(censor_type_i == "clinical"){
        
        pweibull_ab(
            censor_time_i - endpoints_i[1:K + 1],
            censor_time_i - endpoints_i[1:K    ],
            shape = theta$shape_P, scale = theta$scale_P
        )
        
    } else if(censor_type_i == "censored"){
        
        pweibull_ab(
            endpoints_i[1:K    ] - t0, 
            endpoints_i[1:K + 1] - t0, 
            shape = theta$shape_H, scale = theta$scale_H
        )
        
    }
    
    prob_screens <- if(censor_type_i == "screen"){
        (1-theta$beta)^((K-1):0) * theta$beta
    } else if(censor_type_i == "clinical"){
        (1-theta$beta)^((K-1):0)
    } else if(censor_type_i == "censored"){
        (1-theta$beta)^((K-1):0)
    }
    
    prob_tau <- prob_interval * prob_screens
    prob_tau <- prob_tau/sum(prob_tau) # normalize probabilities
    
    return(prob_tau)
    
}

compute_prob_tau <- function(censor_type, censor_time, endpoints, theta, n_cpu, t0){
    mcmapply(
        compute_prob_tau_i,
        censor_type, censor_time, endpoints,
        MoreArgs = list(theta = theta, t0 = t0),
        USE.NAMES = FALSE,
        mc.cores = n_cpu
    )
}

rprop_tau_HP_i <- function(censor_type_i, censor_time_i, endpoints_i, prob_tau_i, theta, t0){
    
    K <- length(prob_tau_i) # number of intervals
    
    # sample interval
    k_new <- sample.int(K, 1, prob = prob_tau_i) # sample interval
    
    # sample tau_HP in chosen interval
    if(censor_type_i == "screen"  ){
        
        sojourn_H_new <- rweibull_trunc(
            endpoints_i[k_new] - t0, endpoints_i[k_new+1] - t0,
            theta$shape_H, theta$scale_H
        )
        tau_HP_new <- t0 + sojourn_H_new
        
    } else if(censor_type_i == "clinical"){
        
        sojourn_P_new <- rweibull_trunc(
            censor_time_i - endpoints_i[k_new+1], censor_time_i - endpoints_i[k_new],
            theta$shape_P, theta$scale_P
        )
        tau_HP_new <- censor_time_i - sojourn_P_new
        
    } else if(censor_type_i == "censored"){
        
        if(k_new == K){ # tau_HP_new > time_censor
            
            tau_HP_new <- Inf
            
        }else{ # tau_HP_new < time_censor
            
            sojourn_H_new <- rweibull_trunc(
                endpoints_i[k_new] - t0, endpoints_i[k_new+1] - t0,
                theta$shape_H, theta$scale_H
            )
            tau_HP_new <- t0 + sojourn_H_new
            
        }
        
    }
    
    return(tau_HP_new)
    
}

rprop_tau_HP <- function(censor_type, censor_time, endpoints, prob_tau, theta, n_cpu, t0){
    mcmapply(
        rprop_tau_HP_i,
        censor_type, censor_time, endpoints, prob_tau,
        MoreArgs = list(theta = theta, t0 = t0),
        USE.NAMES = FALSE,
        mc.cores = n_cpu
    )
}

dlog_prop_tau_HP_i <- function(tau_HP_i, censor_type_i, censor_time_i, endpoints_i, prob_tau_i, theta, t0){
    
    K      <- length(prob_tau_i) # number of intervals
    
    # contribution of k_new
    k_new  <- sum(endpoints_i < tau_HP_i)
    dlog_k <- log(prob_tau_i[k_new])
    
    # contribution of tau_HP_i
    dlog_tau <- if(censor_type_i == "screen"){
        
        sojourn_H <- tau_HP_i - t0
        dweibull_trunc(
            sojourn_H, endpoints_i[k_new] - t0, endpoints_i[k_new+1] - t0,
            theta$shape_H, theta$scale_H, log = T
        )
        
    } else if(censor_type_i == "clinical"){
        
        sojourn_P <- censor_time_i - tau_HP_i
        dweibull_trunc(
            sojourn_P, censor_time_i - endpoints_i[k_new+1], censor_time_i - endpoints_i[k_new],
            theta$shape_P, theta$scale_P, log = T
        )
        
    } else if(censor_type_i == "censored"){
        
        if(k_new == K){ # tau_HP_i > time_censor
            
            0
            
        }else{ # tau_HP_i < time_censor
            
            sojourn_H <- tau_HP_i - t0
            dweibull_trunc(
                sojourn_H, endpoints_i[k_new] - t0, endpoints_i[k_new+1] - t0,
                theta$shape_H, theta$scale_H, log = T
            )
            
        }
        
    }
    
    return(dlog_k + dlog_tau)
    
}

#
# indolent ####

compute_prob_indolent_i <- function(tau_HP_i, censor_type_i, censor_time_i, theta){
    
    L_0 <- dloglik_PI_i(tau_HP_i, indolent_i = 0, censor_type_i, censor_time_i, theta) %>% exp()
    L_1 <- dloglik_PI_i(tau_HP_i, indolent_i = 1, censor_type_i, censor_time_i, theta) %>% exp()
    
    p_i <- L_1 / (L_0 + L_1)
    
}

compute_prob_indolent <- function(tau_HP, censor_type, censor_time, theta, n_cpu){
    mcmapply(
        compute_prob_indolent_i,
        tau_HP, censor_type, censor_time,
        MoreArgs = list(theta = theta),
        USE.NAMES = FALSE,
        mc.cores = n_cpu
    )
}

rprop_indolent_i <- function(censor_type_i, prob_indolent_i){
    
    indolent_i <- if(censor_type_i == "clinical"){
        0
        #}else if(is.infinite(tau_HP_i)){
        #  NA
    }else{
        rbinom(1,1,prob_indolent_i)
    }
    
    return(indolent_i)
    
}

rprop_indolent <- function(censor_type, prob_indolent, n_cpu){
    mcmapply(
        rprop_indolent_i,
        censor_type, prob_indolent,
        #MoreArgs = list(theta = theta),
        USE.NAMES = FALSE,
        mc.cores = n_cpu
    )
}

dlog_prop_indolent_i <- function(indolent_i, censor_type_i, prob_indolent_i){
    
    dlog <- if(censor_type_i == "clinical"){
        0
        #}else if(is.infinite(tau_HP_i)){
        #  0
    }else if(indolent_i){ # Bernoulli pmf
        log(prob_indolent_i)
    }else if(!indolent_i){
        log(1-prob_indolent_i)
    }
    
    return(dlog)
    
}

dlog_prop_indolent <- function(indolent, censor_type, prob_indolent, n_cpu){
    mcmapply(
        dlog_prop_indolent_i,
        indolent, censor_type, prob_indolent,
        #MoreArgs = list(theta = theta),
        mc.cores = n_cpu
    ) %>% sum()
}


#
# psi ####
rprop_psi <- function(theta_cur, epsilon_psi){
    
    # uniform random walk
    psi_prop <- runif(1, theta_cur$psi-epsilon_psi, theta_cur$psi+epsilon_psi)
    
    # reflection on lower bound 0 and upper bound 1 
    if(psi_prop<0){
        psi_prop <- 0 + (0-psi_prop)
    }else if(psi_prop>1){
        psi_prop <- 1 - (psi_prop-1)
    } 
    
    return(psi_prop)
    
}

# latent ####

rprop_latent_i <- function(censor_type_i, censor_time_i, endpoints_i, prob_tau_i, prob_indolent_i, theta, t0){
    
    # propose tau
    tau_HP_new <- rprop_tau_HP_i(censor_type_i, censor_time_i, endpoints_i, prob_tau_i, theta, t0)
    
    # propose indolent
    indolent_new <- rprop_indolent_i(censor_type_i, prob_indolent_i)
    
    #latent_new <- c("tau_HP" = tau_HP_new, "indolent" = indolent_new)
    latent_new <- c(tau_HP_new, indolent_new)
    
    return(latent_new)
    
}

rprop_latent <- function(
        censor_type, censor_time, endpoints, prob_tau, prob_indolent, theta, n_cpu, t0
){
    mcmapply(
        rprop_latent_i,
        censor_type, censor_time, endpoints, prob_tau, prob_indolent,
        MoreArgs = list(theta = theta, t0 = t0),
        USE.NAMES = FALSE,
        mc.cores = n_cpu
    )
}

dlog_prop_latent_i <- function(
        tau_HP_i, indolent_i, censor_type_i, censor_time_i, endpoints_i,
        prob_tau_i, prob_indolent_i, theta, t0
){
    
    # dlog tau
    dlog_prop_tau <- dlog_prop_tau_HP_i(tau_HP_i, censor_type_i, censor_time_i, endpoints_i, prob_tau_i, theta, t0)
    
    # dlog indolent
    dlog_prop_I <- dlog_prop_indolent_i(indolent_i, censor_type_i, prob_indolent_i)
    
    dlog_prop <- dlog_prop_tau + dlog_prop_I
    return(dlog_prop)
    
}

dlog_prop_latent <- function(
        tau_HP, indolent, censor_type, censor_time, endpoints, prob_tau, prob_indolent, theta, n_cpu, t0
){
    mcmapply(
        dlog_prop_latent_i,
        tau_HP, indolent, censor_type, censor_time, endpoints, prob_tau, prob_indolent,
        MoreArgs = list(theta = theta, t0 = t0),
        USE.NAMES = FALSE,
        mc.cores = n_cpu
    ) %>% sum()
}


#
# M-H tau ####

MH_tau_i <- function(tau_HP_cur_i, indolent_i, censor_type_i, censor_time_i, age_screen_i, endpoints_i, n_screen_positive_i, theta, t0){
    
    # propose new latent
    prob_tau_i          <- compute_prob_tau_i(censor_type_i, censor_time_i, endpoints_i, theta, t0)
    tau_HP_new_i        <- rprop_tau_HP_i(censor_type_i, censor_time_i, endpoints_i, prob_tau_i, theta, t0)
    
    # M-H acceptance ratio
    dlog_prop_cur_i <- dlog_prop_tau_HP_i(tau_HP_cur_i, censor_type_i, censor_time_i, endpoints_i, prob_tau_i, theta, t0)
    dlog_prop_new_i <- dlog_prop_tau_HP_i(tau_HP_new_i, censor_type_i, censor_time_i, endpoints_i, prob_tau_i, theta, t0)
    dlog_lik_cur_i  <- dloglik_tau_i(tau_HP_cur_i, indolent_i, censor_type_i, censor_time_i, age_screen_i, n_screen_positive_i, theta, t0)
    dlog_lik_new_i  <- dloglik_tau_i(tau_HP_new_i, indolent_i, censor_type_i, censor_time_i, age_screen_i, n_screen_positive_i, theta, t0)
    
    MH_logratio_i <- dlog_lik_new_i - dlog_lik_cur_i + dlog_prop_cur_i - dlog_prop_new_i
    
    out <- if(runif(1) < exp(MH_logratio_i)){ # accept new latent data
        c(tau_HP_new_i, TRUE)
    }else{ # keep current latent data
        c(tau_HP_cur_i, FALSE)
    }
    
    return(out)
    
}

MH_tau <- function(tau_HP_cur, indolent, censor_type, censor_time, age_screen, endpoints, n_screen_positive, theta, n_cpu, t0){
    mcmapply(
        MH_tau_i,
        tau_HP_cur, indolent, censor_type, censor_time, age_screen, endpoints, n_screen_positive,
        MoreArgs = list(theta = theta, t0 = t0),
        USE.NAMES = FALSE,
        mc.cores = n_cpu
    )
}

#
# M-H (psi, indolent) ####

rprop_dlog_indolent_i <- function(tau_HP_i, censor_type_i, censor_time_i, theta){
    
    prob_indolent_new_i <- compute_prob_indolent_i(tau_HP_i, censor_type_i, censor_time_i, theta)
    indolent_new_i      <- rprop_indolent_i(censor_type_i, prob_indolent_new_i)
    dlog_prop_new_i     <- dlog_prop_indolent_i(indolent_new_i, censor_type_i, prob_indolent_new_i)
    #dlog_lik_new_i      <- dloglik_PI_i(tau_HP_i, indolent_new_i, censor_type_i, censor_time_i, theta)
    
    out <- c(indolent_new_i, dlog_prop_new_i)#, dlog_lik_new_i)
    return(out)
    
}

rprop_dlog_indolent <- function(tau_HP, censor_type, censor_time, theta, n_cpu){
    
    out_mc <- mcmapply(
        rprop_dlog_indolent_i,
        tau_HP, censor_type, censor_time,
        MoreArgs = list(theta = theta),
        USE.NAMES = FALSE,
        mc.cores = n_cpu
    )
    
    indolent_new           <- out_mc[1,]
    dlog_prop_indolent_new <- out_mc[2,] %>% sum()
    #dlog_lik_new           <- out_mc[3,] %>% sum()
    
    out <- list(
        indolent_new=indolent_new, 
        dlog_prop_indolent_new=dlog_prop_indolent_new) 
        #dlog_lik_new=dlog_lik_new)
    return(out)
    
}



MH_psi_indolent <- function(
        theta_cur, indolent_cur, tau_HP, censor_type, censor_time, epsilon_psi, prior, AFS, n_AFS, t0, n_cpu
){
    
    # propose psi
    psi_new   <- rprop_psi(theta_cur, epsilon_psi) # symmetric proposal
    theta_new <- add_psi(theta_cur, psi_new)
    
    # prior density
    dlog_prior_cur <- dbeta(theta_cur$psi, prior$a_psi, prior$b_psi, log = TRUE)
    dlog_prior_new <- dbeta(theta_new$psi, prior$a_psi, prior$b_psi, log = TRUE)
    
    # propose indolent and compute log proposal density
    out          <- rprop_dlog_indolent(tau_HP, censor_type, censor_time, theta_new, n_cpu)
    indolent_new <- out[["indolent_new"]]
    
    # proposal density
    dlog_prop_indolent_new <- out[["dlog_prop_indolent_new"]]
    prob_indolent_cur      <- compute_prob_indolent(tau_HP, censor_type, censor_time, theta_cur, n_cpu)
    dlog_prop_indolent_cur <- dlog_prop_indolent(indolent_cur, censor_type, prob_indolent_cur, n_cpu)
    
    # log likelihood
    dlog_lik_cur <- dloglik_psi(tau_HP, indolent_cur, censor_type, censor_time, theta_cur, AFS, n_AFS, t0, n_cpu)
    dlog_lik_new <- dloglik_psi(tau_HP, indolent_new, censor_type, censor_time, theta_new, AFS, n_AFS, t0, n_cpu)
    
    # M-H acceptance ratio
    (MH_logratio <- (dlog_lik_new + dlog_prior_new - dlog_prop_indolent_new) - 
        (dlog_lik_cur + dlog_prior_cur - dlog_prop_indolent_cur))
    
    out <- if(runif(1) < exp(MH_logratio)){ # accept new values
        list(indolent = indolent_new, theta = theta_new, accept = TRUE)
    } else { # keep current values
        list(indolent = indolent_cur, theta = theta_cur, accept = FALSE)
    }
    
    return(out)
    
}

#
# MCMC ####

MCMC <- function(
        d_obs_screen, d_obs_censor,
        theta_0, prior, 
        epsilon_rate_H, epsilon_rate_P, epsilon_psi,
        t0,
        M = 100, thin = 1,
        n_cpu = 1, verbose = FALSE
        
){
    
    {
        
        # setup
        n_obs <- nrow(d_obs_censor)
        
        M_thin <- M / thin
        RATE_H <- RATE_P <- BETA <- PSI              <- numeric(M_thin)
        ACCEPT_PSI <- ACCEPT_RATE_H <- ACCEPT_RATE_P <- numeric(M_thin)
        TAU_HP <- INDOLENT <- ACCEPT_LATENT <- matrix(nrow = M_thin, ncol = n_obs)
        
        # pre-process data
        theta                   <- theta_0  
        censor_type             <- d_obs_censor$censor_type
        n_screen_positive       <- censor_type=="screen"
        n_screen_positive_total <- sum(d_obs_screen$screen_detected)
        censor_time             <- d_obs_censor$censor_time
        d_obs_screen_tbl        <- d_obs_screen %>% nest(screens = screen_id:screen_detected)
        screens                 <- d_obs_screen_tbl$screens
        age_screen              <- screens %>% map(~ .[["age_screen"]])
        endpoints               <- compute_endpoints(age_screen, censor_type, censor_time, n_cpu, t0)
        
        d_AFS <- count(d_obs_censor, AFS)
        AFS   <- d_AFS$AFS 
        n_AFS <- d_AFS$n 
    
        # initialisation of latent variables
        prob_tau_0      <- compute_prob_tau(censor_type, censor_time, endpoints, theta, n_cpu, t0)
        tau_HP          <- rprop_tau_HP(censor_type, censor_time, endpoints, prob_tau_0, theta, n_cpu, t0)
        prob_indolent_0 <- compute_prob_indolent(tau_HP, censor_type, censor_time, theta, n_cpu)
        indolent        <- rprop_indolent(censor_type, prob_indolent_0, n_cpu)
        
        # MCMC
        tic()
        if(verbose)  time_start <- Sys.time()
        for(m in 1 : M){
            
            # M-H rate_H
            out_rate_H    <- MH_rate_H(tau_HP, theta, prior, censor_time, t0, AFS, n_AFS, n_cpu, epsilon_rate_H)
            theta         <- out_rate_H[["theta" ]]
            accept_rate_H <- out_rate_H[["accept"]]
            
            # M-H rate_P
            out_rate_P    <- MH_rate_P(tau_HP, theta, prior, indolent, censor_type, censor_time, t0, AFS, n_AFS, n_cpu, epsilon_rate_P)
            theta         <- out_rate_P[["theta" ]]
            accept_rate_P <- out_rate_P[["accept"]]
            
            # Gibbs beta
            theta <- gibbs_beta(tau_HP, theta, screens, prior, n_screen_positive_total)
            
            # update (psi, indolent)
            out_psi_indolent <- MH_psi_indolent(theta, indolent, tau_HP, censor_type, censor_time, epsilon_psi, prior, AFS, n_AFS, t0, n_cpu)
            indolent         <- out_psi_indolent[["indolent"]]
            theta            <- out_psi_indolent[["theta"]]
            accept_psi       <- out_psi_indolent[["accept"]]
            
            # update tau_HP
            out_tau       <- MH_tau(tau_HP, indolent, censor_type, censor_time, age_screen, endpoints, n_screen_positive, theta, n_cpu, t0)
            tau_HP        <- out_tau[1,]
            accept_latent <- out_tau[2,]
            
            # save output
            if(m %% thin == 0){
                
                m_thin <- m %/% thin
                
                RATE_H  [m_thin ] <- theta$rate_H
                RATE_P  [m_thin ] <- theta$rate_P
                PSI     [m_thin ] <- theta$psi
                BETA    [m_thin ] <- theta$beta
                TAU_HP  [m_thin,] <- tau_HP
                INDOLENT[m_thin,] <- indolent
                ACCEPT_LATENT[m_thin,] <- accept_latent
                ACCEPT_PSI   [m_thin ] <- accept_psi
                ACCEPT_RATE_H[m_thin ] <- accept_rate_H
                ACCEPT_RATE_P[m_thin ] <- accept_rate_P
                
            }
            
            # print progress
            if(verbose){
                time_current <- Sys.time()
                runtime <- time_current - time_start
                estimated_remaining_time <- runtime * (M-m)/m
                time_end <- time_current + estimated_remaining_time
                print(paste0(
                    "n=", n_obs, " -- ", 
                    m, "/", M, " = ", round(100*m/M, 2), "% -- ",
                    "estimate end: ", time_end
                    ))
            }
            
        }
        
        tictoc  <- toc()
        runtime <- tictoc$toc - tictoc$tic
    } # end runtime
    
    # output
    THETA <- list(
        RATE_H=RATE_H, RATE_P=RATE_P, BETA=BETA, PSI=PSI
    )
    epsilon <- list(
        epsilon_rate_H=epsilon_rate_H, epsilon_rate_P=epsilon_rate_P, epsilon_psi=epsilon_psi
    )
    ACCEPT <- list(
        ACCEPT_RATE_H=ACCEPT_RATE_H, ACCEPT_RATE_P=ACCEPT_RATE_P, 
        ACCEPT_LATENT=ACCEPT_LATENT, ACCEPT_PSI=ACCEPT_PSI
    )
    
    out <- list(
        THETA=THETA, TAU_HP=TAU_HP, INDOLENT=INDOLENT,
        ACCEPT=ACCEPT, epsilon=epsilon, runtime=runtime
    )
    return(out)
    
}
