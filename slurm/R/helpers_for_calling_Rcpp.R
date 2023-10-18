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

#' Transform age at time of screening data to intervals.
#' 
#' @noRd
#' @param data.obj A list object. List contains element `$screen`, `$censored`,
#'   and `$clinical`. Each of these elements is itself a list containing
#'   `$censor_type`, the censoring type (1, 2, 3), `$n`, the number of cases,
#'   `$censor_time`, the censoring times, `$ages_screen` a list describing
#'   the vectorized ages at time of screening, and `$n_screen_positive`
#'   an indicator of positive screening.
#' @param t0 A scalar numeric object. The study starting time point.
#' 
#' @return A list with elements `$screen`, `$censored`, and `$clinical`. Each
#'   of these elements is itself a list containing `$values`, a vector of 
#'   endpoints; `$starts` a vector of length n giving the first index of 
#'   `$values` pertaining to each case; `$ends`, a vector of length n giving 
#'   the last index of `$values` pertaining to each case; and `$lengths`, a 
#'   vector of length n giving the number of values of `$values` that pertain 
#'   to each individual. Note that the `$start` and `$stop` values are prepared 
#'   for use in C++, which indexes from 0.
#' 
#' @keywords internal
compute_endpoints_cpp <- function(data.obj, t0) {
    
    res <- list()
    
    screen <- lapply(seq_len(data.obj$screen$n),
                     FUN = function(i, data, t0) {
                         idx = {data$ages_screen$starts[i] + 1L}:{data$ages_screen$ends[i] + 1L}
                         c(t0, data$ages_screen$values[idx])
                     },
                     data = data.obj$screen, t0 = t0)
    
    res$screen <- list("values"  = unlist(screen),
                       "starts"  = head(c(1L, cumsum(lengths(screen)) + 1L), -1L) - 1L,
                       "ends"    = cumsum(lengths(screen)) - 1L,
                       "lengths" = lengths(screen))
    
    clinical <- lapply(seq_len(data.obj$clinical$n),
                       FUN = function(i, data, t0) {
                           idx = {data$ages_screen$starts[i] + 1L}:{data$ages_screen$ends[i] + 1L}
                           c(t0, 
                             data$ages_screen$values[idx],
                             data$censor_time[i]) |> unique()
                       },
                       data = data.obj$clinical, t0 = t0)
    
    res$clinical <- list("values"  = unlist(clinical),
                         "starts"  = head(c(1L, cumsum(lengths(clinical)) + 1L), -1L) - 1L,
                         "ends"    = cumsum(lengths(clinical)) - 1L,
                         "lengths" = lengths(clinical))
    
    censored <- lapply(seq_len(data.obj$censored$n),
                       FUN = function(i, data, t0) {
                           idx = {data$ages_screen$starts[i] + 1L}:{data$ages_screen$ends[i] + 1L}
                           c(t0, 
                             data$ages_screen$values[idx],
                             data$censor_time[i], Inf) |> unique()
                       },
                       data = data.obj$censored, t0 = t0)
    
    res$censored <- list("values"  = unlist(censored),
                         "starts"  = head(c(1L, cumsum(lengths(censored)) + 1L), -1L) - 1L,
                         "ends"    = cumsum(lengths(censored)) - 1L,
                         "lengths" = lengths(censored))
    res
}

# MCMC ####

MCMC_cpp <- function(
        d_obs_screen, d_obs_censor,
        theta_0, prior, 
        epsilon_rate_H, epsilon_rate_P, epsilon_psi,
        t0,
        M = 100, thin = 1
        ) {
    
    M_thin <- M / thin
    
    # setup
    n_obs <- nrow(d_obs_censor)
    
    # pre-process data
    theta                   <- theta_0  
    
    n_screen_positive_total <- sum(d_obs_screen$screen_detected)
    
    d_obs_screen_tbl        <- d_obs_screen %>% nest(screens = screen_id:screen_detected)
    screens                 <- d_obs_screen_tbl$screens
    age_screen              <- screens %>% map(~ .[["age_screen"]])
    
    # internal function to vectorize the a jagged list
    #
    # `age.screen` is a list. Each element of the list provides the ages at
    #   the time of screening for a single participant.
    # The returned list contains: `$values`, a numeric vector, the input ages 
    #   vectorized; `$starts`, an integer vector of length n, each element, i,
    #   is the first index of `$values` pertaining to participant i; `$ends`,
    #   an integer vector of length n, each element, i, is the last index of
    #   `$values` pertaining to participant i; and `$lengths`, an integer
    #   vector of length n, each element, i, is the number of elements of
    #   `$value` that pertain to participant i. 
    # NOTE: The indices provided in `$start` and `$stop` are generated for use 
    #   in C++, which indexes vectors starting at 0. To use in R, must add
    #   1L.
    .defineAges <- function(age.screen) {
        list("values"  = unlist(age.screen),
             "starts"  = head(c(1L, cumsum(lengths(age.screen)) + 1L), -1L) - 1L,
             "ends"    = cumsum(lengths(age.screen)) - 1L,
             "lengths" = lengths(age.screen))
    }
    
    screens  <- d_obs_censor$censor_type == "screen"
    censored <- d_obs_censor$censor_type == "censored"
    clinical <- d_obs_censor$censor_type == "clinical"
 
    data.obj <- list("screen" = list("censor_type" = 1L,
                                     "n" = sum(screens),
                                     "censor_time" = d_obs_censor$censor_time[screens],
                                     "ages_screen" = .defineAges(age_screen[screens]),
                                     "n_screen_positive" = rep(1L, sum(screens))), 
                     "censored" = list("censor_type" = 2L,
                                       "n" = sum(censored),
                                       "censor_time" = d_obs_censor$censor_time[censored],
                                       "ages_screen" = .defineAges(age_screen[censored]),
                                       "n_screen_positive" = rep(0L, sum(censored))),
                     "clinical" = list("censor_type" = 3L,
                                       "n" = sum(clinical),
                                       "censor_time" = d_obs_censor$censor_time[clinical],
                                       "ages_screen" = .defineAges(age_screen[clinical]),
                                       "n_screen_positive" = rep(0L, sum(clinical))))
    
    endpoints <- compute_endpoints_cpp(data.obj, t0)
    
    # add endpoint information to data.obj
    data.obj$screen$endpoints <- endpoints$screen
    data.obj$censored$endpoints <- endpoints$censored
    data.obj$clinical$endpoints <- endpoints$clinical
    rm(endpoints)
    
    # identify unique ages at first screening and their frequency within
    # the provided data
    # STH - the counting makes me wonder if we are expecting `age` to be
    #   an integer value. If an integer, are we truncating or rounding?
    d_AFS <- count(d_obs_censor, AFS)
    AFS   <- d_AFS$AFS 
    n_AFS <- d_AFS$n 
    
    # initialization of latent variables -- these functions are defined in C++
    
    prob_tau_0 <- compute_prob_tau_List(data.obj, theta, t0)
    age_at_tau_hp_hats <- rprop_age_at_tau_hp_hat_List(data.obj, prob_tau_0, theta, t0)
    # if(need_idolent)
    prob_indolent_0 <- compute_prob_indolent_List(data.obj, age_at_tau_hp_hats, theta)
    #else
    #prob_indolent_0 = list(1)
    
    # if(need_idolent)
    indolents <- rprop_indolent_List(data.obj, prob_indolent_0)
    #else
    #indolent = rep(0, n) # should be 0 (progressive) not 1 (indolent)
    
    
    tic()
    out = MCMC_cpp_internal(data_objects = data.obj, 
             indolents = indolents, 
             prior = prior, 
             age_at_tau_hp_hats = age_at_tau_hp_hats, 
             theta = theta, 
             epsilon_rate_H = epsilon_rate_H, 
             epsilon_rate_P = epsilon_rate_P, 
             epsilon_psi = epsilon_psi,
             t0 = t0, 
             M = M, 
             thin = thin, 
             M_thin = M_thin, 
             n_obs = n_obs,
             n_screen_positive_total = n_screen_positive_total, 
             AFS = AFS, 
             n_AFS = n_AFS)
    runtime = toc()
    
    out[["runtime"]] = runtime$toc - runtime$tic
    
    # re-order the columns of output, so they match the order of the input
    order_cols = c(which(screens), which(censored), which(clinical))
    out[["age_at_tau_hp_hat"]][,order_cols] = out[["age_at_tau_hp_hat"]]
    out[["INDOLENT"]][,order_cols] = out[["INDOLENT"]]
    out[["ACCEPT"]][["ACCEPT_LATENT"]][,order_cols] = out[["ACCEPT"]][["ACCEPT_LATENT"]]
    
    return(out)
}
