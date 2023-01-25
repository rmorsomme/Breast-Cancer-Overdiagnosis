source("helpers.R")
library(tidyverse)

{
  theta <- list(
    rate_H = 0.0005, shape_H = 2,
    rate_P = 0.0015, shape_P = 4,
    beta   = 0.814 , psi     = 0.75 # psi = 0.045
  )
  theta <- update_scales(theta)
  
  age_start <- 40
  
  ages_screen <- seq(age_start, 150, by = 3) # seq(age_start + rpois(), 150, by = 3)
  n_screen <- length(ages_screen)
  
  d_process <- tibble(
    person_id = numeric(n), age_start = numeric(n), 
    sojourn_H = numeric(n), tau_HP    = numeric(n),
    indolent  = numeric(n),
    sojourn_P = numeric(n), tau_PC    = numeric(n)
  )
  d_obs_screen <- tibble(
    person_id  = numeric(), screen_id       = numeric(), 
    age_screen = numeric(), screen_detected = numeric()
  )
  d_obs_censor <- tibble(
    person_id = numeric(n), censor_type = character(n), censor_time = numeric(n)
  )
}

#
# Biological process ####

for(i in 1 : n){ # person i
  
  # Process
  sojourn_H <- draw_sojourn_H(theta)
  tau_HP    <- age_start + sojourn_H
  indolent  <- rbernoulli(1, theta$psi)
  sojourn_P <- if(!indolent) { draw_sojourn_P(theta) } else { Inf }
  tau_PC    <- tau_HP    + sojourn_P
  
  d_process[i, ] <- tibble(
    i, age_start, sojourn_H, tau_HP, indolent, sojourn_P, tau_PC
  )
}

# Observation process ####

for(i in 1 : n){  print(i) # person i
  
  tau_HP    <- d_process$tau_HP[i]
  tau_PC    <- d_process$tau_PC[i]
  age_death <- min(40 + rexp(1,1/20), 120)
  
  for(j in 1 : length(ages_screen)){ # screen j
    
    age_screen <- ages_screen[j]
    
    if(min(age_death, tau_PC) < age_screen){
      if(tau_PC < age_death){ # clinical cancer first
        d_obs_censor[i, ] <- tibble(i, "clinical", tau_PC)
      }else if(age_death < tau_PC){ # death first
        d_obs_censor[i, ] <- tibble(i, "censored", age_death)
      }
      break
    }
    
    compartment     <- compute_compartment(age_screen, tau_HP)
    screen_detected <- screen_result(compartment, theta)
    
    d_obs_screen    <- d_obs_screen %>% 
      add_row(
        person_id = i, screen_id = j,
        age_screen = age_screen,
        screen_detected = screen_detected
      )
    
    if(screen_detected){ # positive screen first
      d_obs_censor[i, ] <- tibble(i, "screen", age_screen)
      break
    }
    
  }
  
}