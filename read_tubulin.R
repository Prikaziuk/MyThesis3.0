# read and process reference (Tubulin)

library(tidyverse)
library(readxl)
library(reshape2)
source('handy.R')

dataPath <- 'D://AGC review/MyThesis3.0/'

read_tubulin <- function(dataPath) {
  
  tub.ABC_path <- paste(dataPath, 'TUB - A, B, C.xlsx', sep = '/')
  
  EXPERIMENTS = c("12", "24", "72", "12r", "24r","72r")
  TYPES = c("zero", "normoxia", "anoxia", "reaer_1", "reaer_24")
  
  tubulin <- NULL
  
  for (exp in EXPERIMENTS) {
    in_tubulin <- read_xlsx(tub.ABC_path, range = cell_limits(c(2, NA), c(NA, 6)), sheet = exp, 
                            col_names = c('gene_name', TYPES)
    ) %>%
      mutate( 
        experiment = exp,
        replicate = rep(c('A', 'B', 'C'), each = 6),
        gene_name = 'Tub'
      ) %>%
      melt(
        id.vars=c('gene_name', 'replicate', 'experiment'), 
        variable.name='types', 
        value.name = 'Ct'
      ) %>%
      group_by(
        gene_name, experiment, types, replicate
      ) %>%
      nest(.key = 'rawCt') %>%
      mutate (
        rawCt = rawCt %>%
          map(replace.outliers_df),
        Ct_mean = rawCt %>%
          map_dbl(colMeans, na.rm = T),
        sd_Ct = rawCt %>%
          map_dbl(colSds)
      )
    
    tubulin <- rbind(tubulin, in_tubulin)
  }
  
  RDSPath <- paste(dataPath, 'RDS_all/tubulin.RDS', sep ='/')
  saveRDS(tubulin, file = RDSPath)
}
