
#
# Take a closer look at the data

d_process
d_obs_screen
d_obs_censor
count(d_obs_censor, censor_type) # (0.8-0.02-0.18)
count(d_process, indolent)

i <- which(d_obs_censor$censor_type == "clinical")[1]
d_process    %>% filter(person_id == i)
d_obs_censor %>% filter(person_id == i)
d_obs_screen %>% filter(person_id == i)

