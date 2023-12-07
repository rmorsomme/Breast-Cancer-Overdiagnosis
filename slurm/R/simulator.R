#
# Helper functions ####
draw_sojourn_H <- function(theta)
  rweibull(1, shape = theta$shape_H, scale = theta$scale_H)

draw_sojourn_P <- function(theta)
  rweibull(1, shape = theta$shape_P, scale = theta$scale_P)

compute_compartment <- function(t, tau_HP){
  case_when(t < tau_HP ~ "H",
            TRUE       ~ "P")
}

screen_result <- function(compartment, theta) {
  if(compartment == "H")  return(FALSE)
  if(compartment == "P")  return(runif(1) < theta$beta)
}

#
# Setup ####
d_process <- tibble(
    person_id = numeric(n), sojourn_H = numeric(n), tau_HP = numeric(n),
    indolent  = numeric(n), sojourn_P = numeric(n), tau_PC = numeric(n)
    )

d_obs_screen <- tibble(
    person_id  = numeric(), screen_id       = numeric(), 
    age_screen = numeric(), screen_detected = numeric()
    )

d_obs_censor <- tibble(
    person_id = numeric(), censor_type = character(), censor_time = numeric(), AFS = numeric()
    )

#
## Biological process ####
for(i in 1 : n){ # person i
    
    sojourn_H <- draw_sojourn_H(theta)
    tau_HP    <- t0 + sojourn_H
    indolent  <- runif(1) < theta$psi
    sojourn_P <- if(!indolent) { draw_sojourn_P(theta) } else { Inf }
    tau_PC    <- tau_HP    + sojourn_P
  
    d_process[i, ] <- tibble(i, sojourn_H, tau_HP, indolent, sojourn_P, tau_PC)
}

#
## Observation process ####
for(i in 1 : n){  if(i%%1e3==0)  print(paste0(i, "/", n)) # person i
  
    AFS       <- sample(40:80, 1, prob = exp(-(40:80)/5)) # age at first screen
    age_death <- min(AFS + rexp(1, 1/5), 100)
  
    tau_HP    <- d_process$tau_HP[i]
    tau_PC    <- d_process$tau_PC[i]
    
    if(tau_PC < AFS)  next # we do not observe individuals that develop a clinical cancer before their first screen.
    
    age_screen <- AFS
    j <- 1
    
    while(age_screen < 1e4){ # while-loop will stop before age_screen > 150, see the three `break` statements
        
      # do a screen
      compartment     <- compute_compartment(age_screen, tau_HP)
      screen_detected <- screen_result(compartment, theta) # screen outcome
      
      # add the screen to the data
      d_obs_screen    <- d_obs_screen %>% 
        add_row(
          person_id = i, screen_id = j,
          age_screen = age_screen,
          screen_detected = screen_detected
        )
      
      # if screen is positive, break
      if(screen_detected){
          d_obs_censor <- d_obs_censor %>% 
              add_row(person_id = i, censor_type  = "screen", censor_time = age_screen, AFS = AFS)
          break
          }
      
      # age of next screen
      inter_screen_interval = 1 + rpois(1, 0.5) # random waiting time until next screen
      age_screen <- age_screen + inter_screen_interval
      j <- j + 1
      
      # if clinical cancer before next screen, break
      if(tau_PC < age_screen){
          d_obs_censor <- d_obs_censor %>% 
              add_row(person_id = i, censor_type  = "clinical", censor_time = tau_PC, AFS = AFS)
          break
          }
      
      # if death before next screen, break
      if(age_death < age_screen){
          d_obs_censor <- d_obs_censor %>% 
              add_row(person_id = i, censor_type  = "censored", censor_time = age_death, AFS = AFS)
          break
          }
    }
}

# d_obs_screen %>% count(person_id  ) %>% count(n) %>% print() # distribution of total number of screens
# d_obs_censor %>% count(censor_type) %>% print()