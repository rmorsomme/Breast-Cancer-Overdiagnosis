draw_tau_HP_rpact <- function(age0, theta, t0 = 40, t1 = 55, t2 = 65){
  
  lambdas <- c(theta$lambda_1, theta$lambda_2, theta$lambda_3)
  
  if(age0 < t0){
    piecewiseSurvivalTime <- c(t0, t1, t2)
    piecewiseLambda       <- lambdas
  } else if(age0 < t1) {
    piecewiseSurvivalTime <- c(age0, t1, t2)
    piecewiseLambda       <- lambdas
  } else if(age0 < t2) {
    piecewiseSurvivalTime <- c(age0, t2)
    piecewiseLambda       <- lambdas[2:3]
  } else {
    piecewiseSurvivalTime <- c(age0)
    piecewiseLambda       <- lambdas[3]
  }
  
  rpact:::.getPiecewiseExponentialRandomNumbersFast(
    n = 1, piecewiseSurvivalTime, piecewiseLambda
  )
  
}

draw_tau_HP <- function(age0, theta, t0 = 40, t1 = 55, t2 = 65){
  
  x <- max(t0, age0) + rexp(1, theta$lambda_1)
  if(x < t1)  return(x)
  
  x <- max(t1, age0) + rexp(1, theta$lambda_2)
  if(x < t2)  return(x)
  
  x <- max(t2, age0) + rexp(1, theta$lambda_3)
  return(x)
}


dlog_tau_HP <- function(x, age0, theta, t0 = 40, t1 = 55, t2 = 65){
  
  # index sets for each of the three intervals
  indices <- 1 : length(x)
  i1      <- indices[x < t1]
  i2      <- indices[between(x, t1, t2)]
  i3      <- indices[t2 < x]
  
  # Now, compute log density for values in each interval in turn
  
  # i1 - values from 1st interval
  i <- i1
  dlog1 <- 
    dexp( # density on 1st interval
      x[i] - pmax(t0, age0[i]), theta$lambda_1,
      log = TRUE
    )
  
  # i2 - values from 2nd interval
  i <- i2
  dlog2 <- 
    pexp( # CDF of 1st interval
      t1 - pmax(t0, age0[i]), theta$lambda_1, 
      lower.tail = FALSE, log.p = TRUE
    ) +
    dexp( # density on 2nd interval
      x[i] - pmax(t1, age0[i]), theta$lambda_2, 
      log = TRUE
    )
  
  # i3 - values from 3rd interval
  i <- i3
  dlog3 <- 
    pexp( # CDF of 1st interval
      t1 - pmax(t0, age0[i]), theta$lambda_1, 
      lower.tail = FALSE, log.p = TRUE
    ) +
    pexp( # CDF of 2nd interval
      t2 - pmax(t1, age0[i]), theta$lambda_2, 
      lower.tail = FALSE, log.p = TRUE
    ) +
    dexp( # density on 2nd interval
      x[i] - pmax(t2, age0[i]), theta$lambda_3, 
      log = TRUE
    )
  
  dlog <- sum(dlog1) + sum(dlog2) + sum(dlog3)
  
  return(dlog)
}