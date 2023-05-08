rm(list=ls())
library(tidyverse)
load("C:/Users/18582/Box/Overdiagnosis_Estimation/BCSC_Data/BCSC_data_ind_40_to_85_2000_2015_v3.RData")

t0 <- 30

# Observed data: endpoint
d_obs_censor <- dat_40_to_85 %>%
    tibble() %>%
    transmute(
        person_id = WomanID,
        AFS = screen_no.1,
        censor_type = case_when(
            case == 0 & screen == 0 ~ "censored",
            case == 1 & screen == 1 ~ "screen"  ,
            case == 1 & screen == 0 ~ "clinical" 
        ),
        censor_time = eventime
    ) %>%
    arrange(person_id)

# Observed data: screens
d_obs_screen <- dat_40_to_85 %>% 
    tibble() %>%
    mutate(
        AFS = screen_no.1
    ) %>% 
    pivot_longer(
        cols=starts_with("screen_no."),
        names_to = "screen_id",
        names_prefix = "screen_no.",
        names_transform = list(screen_id = as.integer),
        values_to = "age_screen",
        values_drop_na = T
    ) %>%
    mutate(screen_detected = (age_screen == eventime) & (screen == 1)) %>%
    select(person_id = WomanID, AFS, screen_id, age_screen, screen_detected) %>%
    arrange(person_id, age_screen)


save(d_obs_censor, d_obs_screen, file = "slurm/data/processed/BCSC_40_to_85.RDATA")


# EDA
#d_obs_censor
#d_obs_screen

#id_screen <- d_obs_censor %>% filter(censor_type=="screen") %>% pull(person_id)
#d_obs_screen %>% filter(person_id %in% id_screen) %>% print(n=1e2)

count(d_obs_censor, censor_type)
count(d_obs_screen, screen_detected)
