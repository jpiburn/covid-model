

source('/covidmodeldata/2020-06-24/00-PARAMS.R')
source('analysis/00-functions.R')

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
#worker_ips <- setdiff(worker_ips, c('172.22.3.211','172.22.3.204', '172.22.3.203'))

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
    split_id = rep(1:(length(cl) - 3), length = n()),
    split_id = as.numeric(split_id),
    split_id = if_else(State == "48", length(cl) + 1, split_id), # texas
    split_id = if_else(State == "51", length(cl) + 2, split_id), # virginia
    split_id = if_else(State == "13", length(cl) + 3, split_id)  # georgia 
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


###################################################################
# Determine which chains to use here....
# 1. if diag is good, if so, use all chains
# 2. else   Check which diagk has the lowest max Rhat 
#           if min_diagk is ok
#              remove k sample from
#           else break
bad_chains <- 
diagnostic_df %>%
  select(-files, -chain_info) %>%
  unnest(cols=diag) %>% 
  filter(Rhat > 1.02) %>%
  arrange(desc(Rhat))

chains_to_drop <- 
  bad_chains %>% 
  select(starts_with("Rhat_drop_")) %>% 
  apply(1, function(x) which(x == min(x)))

bad_chains <- 
  bad_chains %>%
  mutate(
    drop_chain_num = chains_to_drop
  )

diagnostic_df <- mutate(diagnostic_df, good_files = files)

for (i in 1:nrow(bad_chains)) {
  
  state_i <- pull(bad_chains[i, "State"])
  drop_chain_i <- pull(bad_chains[i, "drop_chain_num"])
  
  good_files <- diagnostic_df$good_files[[which(diagnostic_df$State == state_i)]] 
  gooder_files <- good_files[setdiff(1:length(good_files), drop_chain_i)]
  
  diagnostic_df$good_files[[which(diagnostic_df$State == state_i)]] <-  gooder_files
  
  message(glue::glue('chain {drop_chain_i} was removed for State {state_i}'))
}


# boxplot of chain durations
# chain_df <-
#   diagnostic_df %>%
#   select(-files, -diag) %>%
#   unnest(cols=chain_info)
# plot_run_time(diagnostic_df)

save(diagnostic_df, file=DIAG_DF_LOC)
save(diagnostic_df, file="results/2020-06-24/diagnostic.Rdata")

# closing cluster ---------------------------------------------------------
plan(sequential)


