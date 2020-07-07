

source('/covidmodeldata/2020-06-24/00-PARAMS.R')
source('analysis-reformat/00-functions.R')

library(tidyverse)
library(tidyselect)
library(parallel)
library(furrr)
library(rstan)

FIPS='01 02 04 05 06 08 09 10 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 51 53 54 55 56'

###############################################################
# Calculate diagnostics

# # How to remove a thread
# #bad_chains <- list(c('54','3'))
# if(!is.null(bad_chains)){
# for(i in 1:length(bad_chains)){
#   chain <- 
#   from_file <- sprintf('/data/covid/tmp/%s/%s/samples_grw_%s.csv',DATE,bad_chains[[i]][1],bad_chains[[i]][2])
#   to_file <- sprintf('/data/covid/tmp/%s/%s/BAD_samples_grw_%s.csv_BAD',DATE,bad_chains[[i]][1],bad_chains[[i]][2])
#   file.rename(from_file, to_file)
# }
# }




# create cluster connections ----------------------------------------------
worker_ips <- cades_workers$ip_address # from 00-PARAMS.R
worker_ips <- setdiff(worker_ips, c('172.22.3.211','172.22.3.204', '172.22.3.203'))

# Connect and create a cluster
cl <- makeClusterPSOCK(
  worker_ips,
  user = "root",
  rshopts = c(
    "-o", "StrictHostKeyChecking=no",
    "-o", "IdentitiesOnly=yes",
    "-i", SSH_KEY
  ),
  rscript = c("sudo", "docker", "run", "--net=host", "-v", "/covidmodeldata:/covidmodeldata", 
              "rstan/covid", "Rscript"),
  dryrun = FALSE,
  verbose = TRUE
)

plan(list(tweak(cluster, workers = cl), multiprocess))




# group split for cluster -------------------------------------------------
diagnostic_list <- 
  tibble(State = strsplit(FIPS, split=' ')[[1]]) %>%
  mutate(
    files = future_map(.x=State, ~get_sample_paths(.x, DATA_DIR))
    ) %>%
  filter(
    map_lgl(files, function(x) (length(x)>0))
  ) %>%
  mutate(
    split_id = rep(1:18, length = n())
  ) %>% 
  group_by(split_id) %>%
  group_split()


# running diagnostics -----------------------------------------------------
diagnostic_df <- future_map(
  # cluster of workers
  .x = diagnostic_list, 
  
  .f = ~ {
    
    state_df <- .x
    
    state_df <- state_df %>%
      mutate(
        chain_info = future_map(files, get_sampling_params),
        diag       = future_map(
          .x=files,
          ~read_and_diagnose(.x, warmup=WARMUP/NTHIN,
                             par_select=c(starts_with('phi'),
                                          starts_with('tau'),
                                          starts_with('log_lambda'))
          )
        )
      )
    
  state_df
  }
)

diagnostic_df <- bind_rows(diagnostic_df)




# diagnostic_df <- diagnostic_df %>%
#   mutate(
#     diag = future_map(
#       .x=files,
#       ~read_and_diagnose(.x, warmup=WARMUP/NTHIN,
#                          par_select=c(starts_with('phi'),
#                                       starts_with('tau'),
#                                       starts_with('log_lambda'))
#                          )
#       )
#     )

# timing info as well as other chain parameters
# diagnostic_df <- mutate(diagnostic_df,
#     chain_info = future_map(files, get_sampling_params)
#   )

###################################################################
# Determine which chains to use here....
# 1. if diag is good, if so, use all chains
# 2. else   Check which diagk has the lowest max Rhat 
#           if min_diagk is ok
#              remove k sample from
#           else break

diagnostic_df %>%
  select(-files, -chain_info) %>%
  unnest(cols=diag) %>%
  summary()

diagnostic_df %>%
  select(-files, -chain_info) %>%
  unnest(cols=diag) %>% 
  filter(Rhat > 1.02) %>%
  arrange(desc(Rhat)) %>% View()

diagnostic_df %>%
  select(-files, -chain_info) %>%
  unnest(cols=diag) %>%
  arrange(ess_bulk) %>%
  head()

diagnostic_df %>%
  select(-files, -chain_info) %>%
  unnest(cols=diag) %>% View()
  summary()


chain_df <-
  diagnostic_df %>%
  select(-files, -diag) %>%
  unnest(cols=chain_info)

# boxplot of chain durations
plot_run_time(diagnostic_df)

#chain_df %>% filter(State == "34", chain_id != 8) %>% pull(sample_file)


# REPLACE THIS NEXT LINE!!!!!!!!!!!!
diagnostic_df <- mutate(diagnostic_df, good_files=files)

#diagnostic_df[diagnostic_df$State == "37", "good_files"] <- chain_df %>% filter(State == "37", chain_id != 8) %>% pull(sample_file)
# 
diagnostic_df$good_files[[which(diagnostic_df$State == "28")]] <- chain_df %>% filter(State == "28", chain_id != 8) %>% pull(sample_file)

save(diagnostic_df, file=DIAG_DF_LOC)
save(diagnostic_df, file="results/2020-06-24/diagnostic.Rdata")


stopCluster(cl)


