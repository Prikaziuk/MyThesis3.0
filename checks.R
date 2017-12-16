library(tidyverse)

library(readxl)



c_myc_brain <- c(30.72, 30.34, 30.58, 30.34, 30.50, 30.43)
GAPDH_brain <- c(23.70, 23.56, 23.47, 23.65, 23.69, 23.68)

c_myc_kidney <- c(27.06, 27.03, 27.03, 27.10, 26.99, 26.94)
GAPDH_kidney <- c(22.76, 22.61, 22.62, 22.60, 22.61, 22.76)


dCt_brain <- mean(c_myc_brain) - mean(GAPDH_brain)
sd_dCt_brain <- propagate_sd_diff(c_myc_brain, GAPDH_brain)

dCt_kidney <- mean(c_myc_kidney) - mean(GAPDH_kidney)
sd_dCt_kidney <- propagate_sd_diff(c_myc_kidney, GAPDH_kidney)

ddCt <- dCt_kidney - dCt_brain
sd_ddCt <- sd_dCt_kidney


folds <- 2 ** -ddCt
sd_folds <- propagate_sd_power(2, ddCt, sd_ddCt)

# if you want to represent calibrator without error sample should take its error
# propagation sd_dCt to ddCt
sd_ddCt_complete <- sqrt(sd_dCt_kidney ** 2 + sd_dCt_brain ** 2)

sd_folds <- propagate_sd_power(2, ddCt, sd_ddCt_complete)

se_folds <- sd_folds / sqrt(3) * 1.96

folds + sd_folds
folds
folds - sd_folds

mean(c(folds + sd_folds, folds - sd_folds))

# for ddCt error is the same because cov is not defined
# в твоём случае будет определена, потому что будет новое среднее
# биологических повторностей

propagate_sd_diff <- function(minuend, subtrahend) {
  # minuend - vector, Ct or dCt sample
  # subtrahend - vector, Ct or dCt reference
  sqrt(sd(minuend) ** 2 + 
         sd(subtrahend) ** 2 - 
         2 * cov(minuend, subtrahend))
}


propagate_sd_power <- function(base, power, power_sd) {
  # base - <dbl> - 2.0 or efficiency
  # power - <dbl> - ddCt
  # power_sd - <dbl> - propageted sd of dCt
  #                    sd of biological replicates (mean dCt)
  
  # вообще в формуле ещё есть b = -1, но sd = abs() => забыли про него
  (base ** -power) * log(base) * power_sd
}



### figure 2

fos_0 <- c(22.3, 22.0, 21.5)
actin_0 <- c(22.9, 22.3, 22.4)

fos_mean <- mean(fos_0)
actin_mean <- mean(actin_0)
dCt_mean <- fos_mean - actin_mean

dCt <- fos_0 - actin_0
sd_dCt <- propagate_sd_diff(fos_0, actin_0) # == sd_ddCt

ddCt <- dCt - dCt_mean

bio_dCt_mean <- mean(dCt)
bio_dCt_sd <- sd(dCt)

propagate_sd_power(2, mean(ddCt), bio_dCt_sd)

folds <- 2** -ddCt
mean(folds)
sd(folds)

mean(folds) * log(2) * sd_dCt
