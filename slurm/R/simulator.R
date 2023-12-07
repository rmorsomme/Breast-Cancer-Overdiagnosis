sojourn_H <- stats::rweibull(n, shape = theta$shape_H, scale = theta$scale_H)
tau_HP <- t0 + sojourn_H
indolent <- stats::runif(n) < theta$psi
sojourn_P <- stats::rweibull(n, shape = theta$shape_P, scale = theta$scale_P)
sojourn_P[indolent] <- Inf
tau_PC <- tau_HP + sojourn_P

d_obs_screen <- matrix(0.0, nrow = n*30, ncol = 4)
d_obs_censor <- data.frame(person_id = numeric(n),
                           censor_type = character(n),
                           censor_time = numeric(n),
                           AFS = numeric(n))

#
## Observation process ####
k <- 0L
kc <- 0L
AFS <- sample(40L:80L, n, prob = exp(-(40:80)/5), replace = TRUE)
age_death <- pmin(AFS + stats::rexp(n, 1/5), 100.0)

for (i in 1L:n) { 
  if (i %% 1e3 == 0L)  print(paste0(i, "/", n)) # person i
  
  if(tau_PC[i] < AFS[i])  next # we do not observe individuals that develop a clinical cancer before their first screen.
  
  # generate screens until around age_death
  intervals <- 1.0 + stats::rpois(ceiling(age_death[i]) - AFS[i], 0.5)
  ages_screen <- AFS[i] + cumsum(intervals)
  ages_screen <- c(AFS[i], ages_screen)

  for (j in 1L:length(ages_screen)) {
    
    # if clinical cancer before next screen, break
    if (tau_PC[i] < ages_screen[j]) {
      kc <- kc + 1L
      d_obs_censor$person_id[kc] <- i
      d_obs_censor$censor_type[kc] <- "clinical"
      d_obs_censor$censor_time[kc] <- tau_PC[i]
      d_obs_censor$AFS[kc] <- AFS[i]
      break
    }
    
    # if death before next screen, break
    if (age_death[i] < ages_screen[j]) {
      kc <- kc + 1L
      d_obs_censor$person_id[kc] <- i
      d_obs_censor$censor_type[kc] <- "censored"
      d_obs_censor$censor_time[kc] <- age_death[i]
      d_obs_censor$AFS[kc] <- AFS[i]
      break
    }

    if (ages_screen[j] < tau_HP[i]) {
      screen_detected <- FALSE
    } else {
      screen_detected <- stats::runif(1L) < theta$beta
    }
    
    # add the screen to the data
    k <- k + 1L
    d_obs_screen[k, ] <- c(i, j, ages_screen[j], screen_detected)

    # if screen is positive, break
    if (screen_detected) {
      kc <- kc + 1L
      d_obs_censor$person_id[kc] <- i
      d_obs_censor$censor_type[kc] <- "screen"
      d_obs_censor$censor_time[kc] <- ages_screen[j]
      d_obs_censor$AFS[kc] <- AFS[i]
      break
    }
    
  }
}

d_obs_censor <- d_obs_censor[1L:kc, ]
d_obs_screen <- data.frame(d_obs_screen[1L:k, ])
colnames(d_obs_screen) <- c("person_id", "screen_id", "age_screen", "screen_detected")
