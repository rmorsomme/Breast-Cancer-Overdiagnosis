######################################################
# Analytic overdiagnosis calculations
# Version November 14, 2023
# Adaptation of Marc Ryser's original code
######################################################


library(tidyverse)

source("slurm/R/odx/predictODX.R")

t0=30L
age_min <- 50L
age_max <- 74L
screen_freq <- 2L

#-- Load the data

files <- list.files(pth)
t0_vec <- c(0.000000,
            10.000000,
            11.578947, 
            13.894737,
            15.000000,
            16.210526,
            18.526316,
            20.0,
            2.315789,
            20.842105,
            23.157895,
            25.473684,
            27.789474,
            25.0,
            30.105263,
            30.0,
            32.421053,
            34.736842,
            35.000000,
            37.052632,
            39.368421,
            4.631579,
            40.000000,
            41.684211,
            44.000000,
            5.000000,
            6.947368,
            9.263158)

file_name <- "slurm/output/MCMC/BCSC/M=1e+05-AFS_low=50-AFS_upp=74-shape_H=2-shape_P=2-t0=30.RDATA"

load(file_name)
out_odx = predictODX(theta = theta_mcmc, 
                      t0 = t0, 
                      age.min = age_min, 
                      age.max = age_max, 
                      screen.freq = screen_freq,
                      shape.H = 2.0, shape.P = 2.0)


df <- data.frame(t0 = t0_vec, 
                 mean = sapply(result, function(x) x[1L, 1L]), 
                 lower = sapply(result, function(x) x[1L, 2L]), 
                 upper = sapply(result, function(x) x[1L, 3L]))

gg <- ggplot(df, aes(x = t0, y = mean)) + 
  geom_point() + 
  geom_line() + 
  geom_ribbon(data = df, aes(ymin = lower, ymax = upper), linetype = 2, alpha = 0.1)

print(gg)
