#
# Setup ####

{
    d_process <- tibble(
        person_id = numeric(n),
        sojourn_H = numeric(n), tau_HP = numeric(n),
        indolent  = numeric(n),
        sojourn_P = numeric(n), tau_PC = numeric(n)
        )
    
    d_obs_screen <- tibble(
        person_id  = numeric(), screen_id       = numeric(), 
        age_screen = numeric(), screen_detected = numeric()
        )
    
    d_obs_censor <- tibble(
        person_id = numeric(), censor_type = character(), censor_time = numeric(), AFS = numeric()
        )
}

#
## Biological process ####

for(i in 1 : n){ # person i
    
    sojourn_H <- draw_sojourn_H(theta)
    tau_HP    <- t0 + sojourn_H
    indolent  <- runif(1) < theta$psi
    sojourn_P <- if(!indolent) { draw_sojourn_P(theta) } else { Inf }
    tau_PC    <- tau_HP    + sojourn_P
  
    d_process[i, ] <- tibble(
        i, sojourn_H, tau_HP, indolent, sojourn_P, tau_PC
        )

}

#
## Observation process ####

for(i in 1 : n){  print(i) # person i
  
    AFS       <- sample(40:80, 1, prob = exp(-(40:80)/5)) # age at first screen
    age_death <- min(AFS + rexp(1, 1/7), 100)
  
    tau_HP    <- d_process$tau_HP[i]
    tau_PC    <- d_process$tau_PC[i]
    
    if(tau_PC < AFS)  next # we do not observe individuals that develop a clinical cancer before their first screen.
    
    age_screen <- AFS
    j <- 1
    
    while(age_screen < 1e4){ # while-loop will stop before age_screen > 150, see the three `break` statements
        
      compartment     <- compute_compartment(age_screen, tau_HP)
      screen_detected <- screen_result(compartment, theta)
      
      d_obs_screen    <- d_obs_screen %>% 
        add_row(
          person_id = i, screen_id = j,
          age_screen = age_screen,
          screen_detected = screen_detected
        )
      
      if(screen_detected){  # positive screen first
          d_obs_censor <- d_obs_censor %>% 
              add_row(person_id = i, censor_type  = "screen", censor_time = age_screen, AFS = AFS)
          break
          }
      
      age_screen <- age_screen + 1 + rpois(1, 0.25)
      j <- j + 1
      
      if(tau_PC < age_screen){  # clinical cancer first
          d_obs_censor <- d_obs_censor %>% 
              add_row(person_id = i, censor_type  = "clinical", censor_time = tau_PC, AFS = AFS)
          break
          }
      if(age_death < age_screen){  # death first
          d_obs_censor <- d_obs_censor %>% 
              add_row(person_id = i, censor_type  = "censored", censor_time = age_death, AFS = AFS)
          break
          }
      
    }
  
}

count(d_obs_censor, censor_type) %>% print()
count(d_obs_screen, person_id  ) %>% count(n) %>% print()
beep()