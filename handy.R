colSds <- function(single_col_df) {
  # calculates standard deviation as sqrt(variance)
  require(tidyverse)
  
  # df - df with one column
  single_col_df %>% 
    var(na.rm = T) %>%
    sqrt()
}



replace.outliers_df <- function(single_col_df){
  # replaces outliers with NA
  # oultier is above or beyound quartile 
  num.vect <- t(single_col_df)
  q <- quantile(num.vect, na.rm = T)
  lower_quantile <- lo <- q["25%"] - 1.5*(q["75%"] - q["25%"])
  upper_quantil <- up <- q["75%"] + 1.5*(q["75%"] - q["25%"])
  
  num.vect[num.vect < lo | num.vect > up] <-  NA
  return(t(num.vect))
}