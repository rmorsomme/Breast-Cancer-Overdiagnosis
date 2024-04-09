

predictODX <- function(theta, t0, age.min, age.max, screen.freq, 
                       shape.H = 2.0, shape.P = 1.0) {
  
  stopifnot(
    "`theta` must be a list" = !missing(theta) && is.list(theta),
    "`theta` must contain PSI, BETA, RATE_H, RATE_P" = 
      all(c("PSI", "BETA", "RATE_H", "RATE_P") %in% names(theta)),
    "all elements of `theta` must be of the same dimension" = 
      length(theta$PSI) == length(theta$BETA) &&
      length(theta$PSI) == length(theta$RATE_H) &&
      length(theta$PSI) == length(theta$RATE_P),
    "`t0` must be a non-negative scalar" = !missing(t0) &&
      is.numeric(t0) && is.vector(t0) && length(t0) == 1L && t0 >= 0.0,
    "`age.min` must be a non-negative scalar" = !missing(age.min) &&
      is.numeric(age.min) && is.vector(age.min) && length(age.min) == 1L && 
      age.min >= 0.0,
    "`age.max` must be a non-negative scalar >= age.min" = !missing(age.max) &&
      is.numeric(age.max) && is.vector(age.max) && length(age.max) == 1L && 
      age.max >= age.min,
    "`screen.freq` must be a positive scalar" = !missing(screen.freq) &&
      is.numeric(screen.freq) && is.vector(screen.freq) && 
      length(screen.freq) == 1L && screen.freq > 0.0,
    "`shape.H` must be a positive scalar" = !missing(shape.H) &&
      is.numeric(shape.H) && is.vector(shape.H) && 
      length(shape.H) == 1L && shape.H > 0.0,
    "`shape.P` must be a positive scalar" = !missing(shape.P) &&
      is.numeric(shape.P) && is.vector(shape.P) && 
      length(shape.P) == 1L && shape.P > 0.0
  )

  source("slurm/R/odx/predictODX_helpers.R")
  load(file = "slurm/R/odx/OtherCause_birthcohort_1971.RData")
  
  max_n_screens <- (age.max - age.min) / screen.freq + 1L 
  age_at_screens <- seq(age.min, age.max, by = screen.freq)
  
  M <- length(theta$RATE_H)

  screen_ind <- matrix(0.0, nrow = max_n_screens, ncol = M)
  screen_prog <- matrix(0.0, nrow = max_n_screens, ncol = M)
  screen_total <- matrix(0.0, nrow = max_n_screens, ncol = M)
  
  lambda_seq <- matrix(0.0, M, 2L)
  
  cat(paste(0:10, collapse = "----"), "\n", sep = "")
  tick_increment = max(floor(M/50), 1L)
  
  for (m in seq_len(M)) {
    
    if (m %% tick_increment == 0L) cat("=")
    
    # model parameters
    sample_theta <- lapply(theta, "[", m)
    
    # model parameters
    psi <- sample_theta$PSI
    beta <- sample_theta$BETA
    scale_H <- sample_theta$RATE_H^(-1.0 / shape.H)
    scale_P <- sample_theta$RATE_P^(-1.0 / shape.P)
    
    ## The lambda used
    lambda_seq[m, ] <- c(scale_P, shape.P)
    
    for (tt in (0L:(max_n_screens-1)*2)) {
      
      t.vec <- c(t0, seq(from = age.min, to = (age.min + tt), by = screen.freq))
      
      # number of screens
      ni <- length(t.vec) - 1L
      
      Dj <- D_j(psi, beta, shape.H, scale_H, shape.P, scale_P, t.vec,  ni)
      
      # P(SD)
      screen_total[ni, m] <- {Dj$indolent + Dj$progressive} * Dj$norm
      
      # P(SD, I=1)
      screen_ind[ni, m] <- Dj$indolent * Dj$norm
      
      # P(SD, I=0)
      screen_prog[ni, m] <- Dj$progressive * Dj$norm
      
    }
  } 
  cat('\n')
  
  #-- P(screen-detected cancer)
  SD_stat <- numeric(M)
  #-- P(indolent screen-detected)
  ODX_ind_stat <- numeric(M)
  #-- P(progressive overdiagnosed cancer)
  ODX_prog_stat <- numeric(M)
  
  cat(paste(0:10, collapse = "----"), "\n", sep = "")
  
  for (z in seq_len(M)) {
    
    if (z %% tick_increment == 0L) cat("=")
    
    #--Probability to have a screen-detected cancer
    # probability to survive past each respective screening round
    hh <- sapply(age_at_screens, surv.fun, age.min = age.min, 
                 lambda.other = OChaz$Estimate, lambda.ages = OChaz$Age)
    
    # probability of a screen-detected cancer while stimax_n_screens alive
    SD_stat[z] <- sum(screen_total[, z] * hh)
    
    ## Probability to have an indolent screen-detected cancer
    
    # probability of an indolent screen-detected cancer while stimax_n_screens alive
    ODX_ind_stat[z] <- sum(screen_ind[, z] * hh)
    
    ## Probability to have a progressive and overdiagnosed screen-detected cancer
    hh <- sapply(age_at_screens, Prog.odx.program, age.min = age.min, 
                 P.lambda = lambda_seq[z,], 
                 lambda.other = OChaz$Estimate, lambda.ages = OChaz$Age)
    
    ODX_prog_stat[z] <- sum(screen_prog[,z] * hh)
    
  }
  cat('\n')
  
  wrap(ODX_ind_stat, ODX_prog_stat, SD_stat)
}
