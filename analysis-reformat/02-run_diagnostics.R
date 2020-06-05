#ROOT <- '/home/nnagle/Dropbox/students/Piburn/covid-model'
#source('analysis-reformat/00-PARAMS.R')
DATE <- '2020-06-02'
source(file.path('/data/covid/tmp', DATE, '00-PARAMS.R'))
source('analysis-reformat/00-functions.R')

# This uses the future library and furrr package for multithreading.
# I use availableCores()-4 on my own machine
# On my machinge with 32 cores (28 used) it takes about 15 minutes (mostly reading)
#
# The output is a small diagnostic.Rdata file in DATA_DIR

library(tidyverse)
library(tidyselect)
#library(vroom)
library(parallel)
#library(tidyr)
library(furrr)
library(rstan) # The Rhat function is used by one of the readers.

FIPS='01 02 04 05 06 08 09 10 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 51 53 54 55 56'

#DATA_DIR <- '/data/covid/tmp/2020-05-18'
#FIPS = paste(basename(list.dirs(DATA_DIR))[-1], collapse=' ' )

#cl <- makeClusterPSOCK(availableCores() - 70)
cl <- makeClusterPSOCK(NNODES)
plan(cluster, workers = cl)


###############################################################
# Calculate diagnostics

# How to remove a thread
#bad_chains <- list(c('54','3'))
if(!is.null(bad_chains)){
for(i in 1:length(bad_chains)){
  chain <- 
  from_file <- sprintf('/data/covid/tmp/%s/%s/samples_grw_%s.csv',DATE,bad_chains[[i]][1],bad_chains[[i]][2])
  to_file <- sprintf('/data/covid/tmp/%s/%s/BAD_samples_grw_%s.csv_BAD',DATE,bad_chains[[i]][1],bad_chains[[i]][2])
  file.rename(from_file, to_file)
}
}
  
# Evaluate Rhat for the parameter bo_raw
# using all chains
diagnostic_df <- 
  tibble(
    State = strsplit(FIPS, split=' ')[[1]]) %>%
  mutate(
    files = future_map(.x=State, 
                       ~get_sample_paths(.x, DATA_DIR)))
diagnostic_df <- diagnostic_df %>%
  filter(map_lgl(diagnostic_df$files, function(x) (length(x)>0)))

diagnostic_df <- diagnostic_df %>%
  mutate(
    diag = future_map(.x=files, 
                      ~read_and_diagnose(.x, 
                                         warmup=WARMUP/NTHIN, 
                                         par_select=c(starts_with('phi'),
                                                      starts_with('tau'),
                                                      starts_with('log_lambda')))))

# timing info as well as other chain parameters
diagnostic_df <- mutate(diagnostic_df,
    chain_info = future_map(files, get_sampling_params)
  )

###################################################################
# Determine which chains to use here....
# 1. if diag is good, if so, use all chains
# 2. else   Check which diagk has the lowest max Rhat 
#           if min_diagk is ok
#              remove k sample from
#           else break

diagnostic_df %>%
  select(-files) %>%
  unnest(cols=diag) %>%
  summary()

diagnostic_df %>%
  select(-files) %>%
  unnest(cols=diag) %>% 
  filter(Rhat > 1.02) %>%
  arrange(desc(Rhat))

diagnostic_df %>%
  select(-files) %>%
  unnest(cols=diag) %>%
  arrange(ess_bulk) %>%
  head()

# boxplot of chain durations
plot_run_time(diagnostic_df)


# REPLACE THIS NEXT LINE!!!!!!!!!!!!
diagnostic_df <- mutate(diagnostic_df, good_files=files)


save(diagnostic_df, file=DIAG_DF_LOC)

stopCluster(cl)


if (Sys.info()['sysname']=='Windows'){
  system2('taskkill',args=c('/f','/im "Rscript.exe"'), stdout=FALSE)
}
