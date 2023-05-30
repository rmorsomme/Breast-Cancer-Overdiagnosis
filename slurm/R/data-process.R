
rm(list=ls())
library(tidyverse)
library(lubridate)



#
# BCSC data ####
load("C:/Users/18582/Box/Overdiagnosis_Estimation/BCSC_Data/BCSC_data_ind_40_to_85_2000_2015_v3.RData")

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
# TODO: remove AFS from d_obs_screen
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




#
# Swiss data ####

# The goal is to create two data sets
# 1. d_obs_censor with the variables: person_id, AFS (age at first screen, in years), censor_type, censor_time (in years)
# 2. d_obs_screen with the variables: person_id, screen_id, age_screen (in years), screen_detected (binary)
# for an example, see the code above for the BCSC data

## Read in raw data ####
path <- "C:/Users/18582/Box/Overdiagnosis_Estimation/Swiss_data/Pull_3/" # most recent data

d_comp     <- read_delim(file = paste0(path, "compliance.csv"), delim = ";")
d_exam_tmp <- read_delim(file = paste0(path, "mammogram.csv" ), delim = ";") # temporary data set; will be merged with d_comp shortly
d_woman    <- read_delim(file = paste0(path, "woman.csv"     ), delim = ";") %>%
    mutate(
        birthdate      = dmy(birthdate),
        endpointdt     = if_else(endpointdt != "99999999", endpointdt, NA_character_) %>% dmy,
        age_endpoint   = ((endpointdt - birthdate) / 365) %>% as.numeric()
    ) %>%
    select(person_id = nodoss, birthdate, age_endpoint, detec_amended_recoding, endpointdt, endpoint_year)


## Merge data sets ####
all(d_comp$nodoss == d_exam_tmp$nodoss & d_comp$examdt == d_exam_tmp$examdt) # check coherence across d_comp and d_exam_tmp

d_exam <- left_join(
    d_comp     %>% mutate(row_id = 1:nrow(.)),
    d_exam_tmp %>% mutate(row_id = 1:nrow(.)),
    by = c("row_id", "nodoss", "examdt") # have to join by row_id because there are multiple rows with the same nodoss and with examdt=9999999
    ) %>%
    mutate(examdt = if_else(examdt != "99999999", examdt, NA_character_) %>% dmy) %>% 
    filter(!is.na(examdt) & concl != 9) %>% # remove screens with a missing date or a missing BI-RADS
    select(person_id = nodoss, examdt)

d_full <- inner_join(d_exam, d_woman, by = "person_id") %>% # use an inner join because the two data sets are not perfectly sync'ed
    mutate(
        age_screen = ((examdt - birthdate) / 365) %>% as.numeric(),
        censor_type = case_when( # define the type of censoring based on the mode of detection, not the outcome of the screens
            detec_amended_recoding == 0        ~ "clinical",
            detec_amended_recoding %in% c(1,4) ~ "screen",
            detec_amended_recoding %in% c(2,5) ~ "censored", # we will manually censor the individuals that are detected out of program
            detec_amended_recoding == 3        ~ NA_character_
        )
    ) %>%
    filter(
        !is.na(censor_type),
        age_screen <= age_endpoint
        ) %>%
    group_by(person_id) %>%
    mutate(
        AFS             = min(age_screen), # age at first screen
        screen_id       = 1 : n(),
        exam_year_first = examdt %>% min %>% year # calendar year of first exam (will be used to take a subset of the data)
    ) %>%
    ungroup() %>%
    select(person_id, screen_id, age_screen, exam_year_first, age_endpoint, endpoint_year, AFS, censor_type)

count(d_full, exam_year_first) %>% print(n=1e2)

## Take a subset ####
d_full_sub <- d_full %>%
    filter(
        exam_year_first %>% between(2010, 2015), # only keep participants that entered the program between 2010 and 2015
        ! (censor_type == "clinical" & age_screen > age_endpoint - 0.25) # remove screens that are very close (3 months) to the endpoint date for participants with clinical cancer.
        ) %>%
    select(person_id, screen_id, age_screen, age_endpoint, AFS, censor_type) %>%
    group_by(person_id) %>%
    mutate( # generate the screen outcome based on the 
        screen_detected = if_else(
            censor_type == "screen", 
            rep(c(0, 1), c(n()-1, 1)), # (n-1) negative screens and 1 positive screens
            rep(  0    ,   n()      )  # n negative screens
            )
    )%>%
    ungroup()


# Output -- 2 dataframes

# 1. d_obs_censor
# with the variables: person_id, AFS (age at first screen, in years), censor_type, censor_time (in years)
d_obs_censor <- d_full_sub %>%
    group_by(person_id) %>%
    summarize( # arguably, there should be a much more elegant way to wrangle the data which does not make us do this stupid summarize() call; should probably nest the screen and their outcomes in a variable (column) of dataframes
        ALS = max(age_screen), # age at last screen
        AFS = unique(AFS),
        censor_type = unique(censor_type),
        age_endpoint = unique(age_endpoint)
        ) %>%
    mutate(
        censor_time = case_when(
            censor_type == "clinical" ~ age_endpoint,
            censor_type == "screen"   ~ ALS, # or should I use age_endpoint for screen-detected individuals? These two variables do not necessarily agree!
            censor_type == "censored" ~ pmin(ALS + 1.5, age_endpoint)
        )
    ) %>%
    select(-ALS, - age_endpoint)

# 2. d_obs_screen with the variables: person_id, screen_id, age_screen (in years), screen_detected (binary)
d_obs_screen <- d_full_sub %>% select(person_id, screen_id, age_screen, screen_detected)

save(d_obs_censor, d_obs_screen, file = "slurm/data/processed/swiss.RDATA")


count(d_obs_censor, censor_type)
count(d_obs_screen, screen_detected)
