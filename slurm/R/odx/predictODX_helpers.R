### Calculating over-diagnosis contributions

#' Integral f(t; shape.H, scale.H) S(Y-t; shape.P, scale.P) from L:U
#' 
#' @noRd
#' @param U A scalar numeric. The upper integration boundary.
#' @param L A scalar numeric. The lower integration boundary.
#' @param last A scalar numeric. The final screening age.
#' @param shape.H A scalar numeric. The shape parameter for the Weibull of the
#'   healthy compartment.
#' @param scale.H A scalar numeric. The scale parameter for the Weibull of the
#'   healthy compartment.
#' @param shape.P A scalar numeric. The shape parameter for the Weibull of the
#'   pre-clinical compartment.
#' @param scale.P A scalar numeric. The scale parameter for the Weibull of the
#'   pre-clinical compartment.
#'
#' @return A scalar numeric.
#'
#' @keywords internal
.integral <- function(U, L, last, shape.H, scale.H, shape.P, scale.P) {

  integrand <- function(t, last, shape.H, shape.P, scale.H, scale.P){
    stats::dweibull(t, shape.H, scale.H) * 
      stats::pweibull(last - t, shape.P, scale.P, lower.tail = FALSE)
  }
  
  stats::integrate(integrand, lower = L, upper = U,
                   last = last, 
                   shape.H = shape.H, shape.P = shape.P, 
                   scale.H = scale.H, scale.P = scale.P, 
                   stop.on.error = TRUE)$value
  
}

#' Contribution to the normalization for non-indolent
#' 
#' @noRd
#' @param U A scalar numeric. The upper integration boundary.
#' @param L A scalar numeric. The lower integration boundary.
#' @param last A scalar numeric. The final screening age.
#' @param shape.H A scalar numeric. The shape parameter for the Weibull of the
#'   healthy compartment.
#' @param scale.H A scalar numeric. The scale parameter for the Weibull of the
#'   healthy compartment.
#' @param shape.P A scalar numeric. The shape parameter for the Weibull of the
#'   pre-clinical compartment.
#' @param scale.P A scalar numeric. The scale parameter for the Weibull of the
#'   pre-clinical compartment.
#'
#' @return A scalar numeric.
#'
#' @keywords internal
.normalization_not_indolent <- function(U, L,  
                                        shape.H, scale.H, 
                                        shape.P, scale.P) {
  
  stats::pweibull(U, shape.H, scale.H, lower.tail = FALSE) + 
    .integral(U = U, L = L, last = U, 
              shape.H = shape.H, scale.H = scale.H, 
              shape.P = shape.P, scale.P = scale.P)
}

#' Normalization
#' 
#' Accounts for the time lapse between risk onset and first screen.
#' 
#' @noRd
#' @param psi A scalar numeric. The probability of being indolent.
#' @param U A scalar numeric. The upper integration boundary.
#' @param L A scalar numeric. The lower integration boundary.
#' @param last A scalar numeric. The final screening age.
#' @param shape.H A scalar numeric. The shape parameter for the Weibull of the
#'   healthy compartment.
#' @param scale.H A scalar numeric. The scale parameter for the Weibull of the
#'   healthy compartment.
#' @param shape.P A scalar numeric. The shape parameter for the Weibull of the
#'   pre-clinical compartment.
#' @param scale.P A scalar numeric. The scale parameter for the Weibull of the
#'   pre-clinical compartment.
#'
#' @return A scalar numeric.
#'
#' @keywords internal
.normalization <- function(psi, U, L,
                           shape.H, scale.H, shape.P, scale.P) {
  
  psi + (1.0 - psi) * 
    .normalization_not_indolent(U = U, L = L, 
                                shape.H = shape.H, scale.H = scale.H, 
                                shape.P = shape.P, scale.P = scale.P)
}

#' ALL SCREEN-DETECTED CASES likelihood contribution
#'
#' @noRd
#' @param psi A scalar numeric. The probability of indolence.
#' @param shape.H A scalar numeric. The shape parameter for the Weibull of the
#'   healthy compartment.
#' @param scale.H A scalar numeric. The scale parameter for the Weibull of the
#'   healthy compartment.
#' @param shape.P A scalar numeric. The shape parameter for the Weibull of the
#'   pre-clinical compartment.
#' @param scale.P A scalar numeric. The scale parameter for the Weibull of the
#'   pre-clinical compartment.
#' @param time_points A numeric vector object. The age at risk onset and
#'   ages at each screen
#'   (age_risk_onset, age_screen_1, ..., age_screen_t)
#' @param ni An integer object. The number of screens.
#' 
#' @returns A scalar numeric. The contribution to the likelihood for screen
#'   detected cases
#' @keywords internal
D_j <- function(psi, beta, shape.H, scale.H, shape.P, scale.P, time_points, ni) {
    
  # the number of time points 
  l <- ni + 1L
  
  t0 <- time_points[1L]
  AFS <- time_points[2L]
  
  # initializations
  A <- 0.0
  B <- 0.0

  ## sum over previous screens 
  for (k in 2L:l) {
    tmp <- beta * {(1.0 - beta)^{l - k}}
    A <- A +  psi * tmp * 
      {stats::pweibull(time_points[k] - t0, shape.H, scale.H) - 
          stats::pweibull(time_points[k - 1L] - t0, shape.H, scale.H)}
   
    B <- B + (1.0 - psi) * tmp * 
      .integral(U = time_points[k] - t0, L = time_points[k - 1L] - t0,
                last = time_points[l] - t0,
                shape.H = shape.H, scale.H = scale.H,
                shape.P = shape.P, scale.P = scale.P)
  }
  
  norm <- .normalization(psi = psi, U = AFS - t0, L = 0.0,
                         shape.H = shape.H, scale.H = scale.H, 
                         shape.P = shape.P, scale.P = scale.P)
  
  list("indolent" = A, "progressive" = B, "norm" = 1.0 / norm)
}  

###############################################
### Other cause death functions  ##############
###############################################


#' other cause mortality survival function 
#' 
#' @noRd
#' @param age.min A scalar numeric. The age of first screening. 
#'   Note that we are assuming that it is an integer value.
#' @param t A scalar numeric. The time to evaluate the survival function.
#' @param lambda.other A numeric vector. The other cause hazard
#' 
#' @returns A numeric vector.
#' @keywords internal
surv.fun <- function(age.min, t, lambda.other, lambda.ages) {

  if (t < age.min) return(1.0)
  
  bounds <- findInterval(c(age.min, t), lambda.ages)
  lambdas <- lambda.other[bounds[1L]:bounds[2L]]
  deltat <- lambda.ages[bounds[1L]:bounds[2L]]
  deltat <- c(deltat[-1L], t) - deltat
  hh <- drop(crossprod(lambdas, deltat))

  exp(-hh)
}

surv.fun <- Vectorize(surv.fun, vectorize.args = "t")

########### Prog.odx.programmatic

#' @noRd
#' @param age.min A numeric scalar. Age at first screening.
#' @param screen.age A numeric scalar. The screen age.
#' @param lambda A numeric scalar. The modeled risk
#' @param lambda.other A numeric vector. The other cause hazard
#'
#' @returns 
#' @keywords internal
Prog.odx.program <- function(age.min, screen.age, P.lambda, lambda.other, lambda.ages) {
  
  help.f <- function(t) {
    surv.fun(age.min = age.min, t = screen.age + t, lambda.other = lambda.other, 
             lambda.ages = lambda.ages) * 
      stats::dweibull(x = t, shape = P.lambda[2L], scale = P.lambda[1L])
  }
  
  out <- integrate(help.f, lower = 0.0, upper = 500.0)

  surv.fun(age.min, screen.age, lambda.other, lambda.ages) - out$value
}

Prog.odx.program <- Vectorize(Prog.odx.program, vectorize.args = "screen.age")



########### Wrapper for the summary statistics

wrap <- function(ODX.ind.stat, ODX.prog.stat, SD.stat){
  
  
  
  #overall
  a<-c(round(100.0 * mean((ODX.ind.stat + ODX.prog.stat) / SD.stat), 1),
       round(100.0 * stats::quantile((ODX.ind.stat + ODX.prog.stat) / SD.stat, 
                                     probs = c(0.025, 0.5, 0.975)), 1))
  #indolent
  b<-c(round(100.0 * mean(ODX.ind.stat / SD.stat), 1),
       round(100.0 * stats::quantile(ODX.ind.stat / SD.stat, 
                                     probs = c(0.025, 0.5, 0.975)), 1))
  #mortality
  c<-c(round(100.0 * mean(ODX.prog.stat / SD.stat), 1),
       round(100.0 * stats::quantile(ODX.prog.stat / SD.stat, 
                                     probs = c(0.025, 0.5, 0.975)), 1))
  
  
  tab <- rbind(a, b, c)
  colnames(tab)[1L] <- "mean"
  tab <- tab[, c(1L, 2L, 4L, 3L)]
  rownames(tab) <- c("total", "indolent", "mortality")
  
  tab
}


