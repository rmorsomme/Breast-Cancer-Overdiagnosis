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

MH_tau_R <- function(tau_HP, indolent, censor_type, censor_time, 
                     age_screen, endpoints, n_screen_positive, theta, t0) {
    
    func <- function(i, tau_HP, indolent, censor_type, censor_time, 
                     age_screen, endpoints, n_screen_positive, theta, t0) {
        MH_tau_i(tau_HP[i], indolent[i], censor_type[i], censor_time[i], age_screen[[i]], 
                 endpoints[[i]], n_screen_positive[i], theta, t0)
    }
    
    result <- future.apply::future_lapply(X = seq_along(tau_HP),
                                          FUN = func,
                                          tau_HP = tau_HP, 
                                          indolent = indolent, 
                                          censor_type = censor_type, 
                                          censor_time = censor_time, 
                                          age_screen = age_screen, 
                                          endpoints = endpoints, 
                                          n_screen_positive = n_screen_positive,
                                          theta = theta, t0 = t0,
                                          future.seed = TRUE)
    
    result_tau <- lapply(result, "[[", 1) |> unlist()
    result_accept <- lapply(result, "[[", 2) |> unlist()
    
    return(list("tau_HP" = result_tau, 
                "accept" = result_accept))
}

#' @noRd
#' @param data.obj A list object. List contains element `$screen`, `$censored`,
#'   and `$clinical`. Each of these elements is itself a list containing
#'   `$endpoints`, `$tau.HP`, `$censor_time`, `$ages_screen`, and `$not.indolent`.
#' @param t0 A scalar numeric object.
#' @return A list containing `$val`, a vector of endpoints, `$starts` a vector
#'   of length n giving the first index pertaining to each case, and `$ends`,
#'   a vector of length n giving the last index pertaining to each case
#' 
#' @keywords internal
compute_endpoints <- function(data.obj, t0) {
    
    res <- list()
    
    screen <- lapply(seq_len(data.obj$screen$n),
                     FUN = function(i, data, t0) {
                         idx = {data$ages_screen$starts[i] + 1}:{data$ages_screen$ends[i] + 1}
                         c(t0, data$ages_screen$values[idx])
                     },
                     data = data.obj$screen, t0 = t0)
    
    res$screen <- list("values"  = unlist(screen),
                       "starts"  = head(c(1L, cumsum(lengths(screen)) + 1L), -1L) - 1L,
                       "ends"    = cumsum(lengths(screen)) - 1L,
                       "lengths" = lengths(screen))
    
    clinical <- lapply(seq_len(data.obj$clinical$n),
                       FUN = function(i, data, t0) {
                           idx = {data$ages_screen$starts[i] + 1}:{data$ages_screen$ends[i] + 1}
                           c(t0, 
                             data$ages_screen$values[idx],
                             data$censor_time[i])
                       },
                       data = data.obj$clinical, t0 = t0)
    
    res$clinical <- list("values"     = unlist(clinical),
                         "starts"  = head(c(1L, cumsum(lengths(clinical)) + 1L), -1L) - 1L,
                         "ends"    = cumsum(lengths(clinical)) - 1L,
                         "lengths" = lengths(clinical))
    
    censored <- lapply(seq_len(data.obj$censored$n),
                       FUN = function(i, data, t0) {
                           idx = {data$ages_screen$starts[i] + 1}:{data$ages_screen$ends[i] + 1}
                           c(t0, 
                             data$ages_screen$values[idx],
                             data$censor_time[i], Inf)
                       },
                       data = data.obj$censored, t0 = t0)
    
    res$censored <- list("values"     = unlist(censored),
                         "starts"  = head(c(1L, cumsum(lengths(censored)) + 1L), -1L) - 1L,
                         "ends"    = cumsum(lengths(censored)) - 1L,
                         "lengths" = lengths(censored))
    res
}

# MCMC ####

MCMC <- function(
        d_obs_screen, d_obs_censor,
        theta_0, prior, 
        epsilon_rate_H, epsilon_rate_P, epsilon_psi,
        t0,
        M = 100, thin = 1,
        n_cpu = NULL, verbose = FALSE) {
    
    M_thin <- M / thin
    
    # setup
    n_obs <- nrow(d_obs_censor)
    
    # pre-process data
    theta                   <- theta_0  
    
    n_screen_positive_total <- sum(d_obs_screen$screen_detected)
    
    d_obs_screen_tbl        <- d_obs_screen %>% nest(screens = screen_id:screen_detected)
    screens                 <- d_obs_screen_tbl$screens
    age_screen              <- screens %>% map(~ .[["age_screen"]])
    
    .defineAges <- function(age_screen) {
        list("values"  = unlist(age_screen),
             "starts"  = head(c(1L, cumsum(lengths(age_screen)) + 1L), -1L) - 1L,
             "ends"    = cumsum(lengths(age_screen)) - 1L,
             "lengths" = lengths(age_screen))
    }
    # screen cases
    screens <- d_obs_censor$censor_type == "screen"
    screen_object = list("censor_type" = 1L,
                         "n" = sum(screens),
                         "censor_time" = d_obs_censor$censor_time[screens],
                         "ages_screen" = .defineAges(age_screen[screens]),
                         "n_screen_positive" = rep(1L, sum(screens)))
    
    # censored cases
    censored <- d_obs_censor$censor_type == "censored"
    censored_object = list("censor_type" = 2L,
                           "n" = sum(censored),
                           "censor_time" = d_obs_censor$censor_time[censored],
                           "ages_screen" = .defineAges(age_screen[censored]),
                           "n_screen_positive" = rep(0L, sum(censored)))
    
    # censored cases
    clinical <- d_obs_censor$censor_type == "clinical"
    
    clinical_object = list("censor_type" = 3L,
                           "n" = sum(clinical),
                           "censor_time" = d_obs_censor$censor_time[clinical],
                           "ages_screen" = .defineAges(age_screen[clinical]),
                           "n_screen_positive" = rep(0L, sum(clinical)))
    
    data.obj <- list("screen" = screen_object, 
                     "censored" = censored_object,
                     "clinical" = clinical_object)
    
    endpoints <- compute_endpoints(data.obj, t0)
    data.obj$screen$endpoints <- endpoints$screen
    data.obj$censored$endpoints <- endpoints$censored
    data.obj$clinical$endpoints <- endpoints$clinical
    
    d_AFS <- count(d_obs_censor, AFS)
    AFS   <- d_AFS$AFS 
    n_AFS <- d_AFS$n 
    
    # initialisation of latent variables
    prob_tau_0      <- compute_prob_tau_List(data.obj, theta, t0)
    tau_HPs          <- rprop_tau_HP_List(data.obj, prob_tau_0, theta, t0)
    prob_indolent_0 <- compute_prob_indolent_List(data.obj, tau_HPs, theta)
    indolents        <- rprop_indolent_List(data.obj, prob_indolent_0)
    
    pa <- proc.time()
    out <- MCMC_cpp(data.obj, indolents, prior, tau_HPs, theta, 
                    epsilon_rate_H, epsilon_rate_P, epsilon_psi,
                    t0, M, thin, M_thin, n_obs,
                    n_screen_positive_total, AFS, n_AFS)
    pb <- proc.time()
    print(pb - pa)
    out
}
