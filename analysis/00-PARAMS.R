DATE <- '2020-05-03' # Date of run
NYT_FILE <- NULL
ACS_FILE <- NULL
DATE_0 <- '2020-03-01' # First date to use
##ROOT = ""
SAMPLES_ROOT <- ifelse(Sys.info()[["nodename"]] == 'quetelet',
                       '/data/covid/tmp',
                       'tmp') # All samples will be stored under {SAMPLES_ROOT}/{DATE}/{STATE}
NITER = 3000
NTHIN = 2
NCHAINS = 4
WARMUP = (NITER/NTHIN)/2
WARMUP = 500 # WARMUP is number of retained (after thinning) samples
NNODES = ifelse(Sys.info()[["nodename"]] == 'quetelet',
                30,
                80) # Number of nodes in cluster
DATA_DIR <- file.path(SAMPLES_ROOT,DATE)
DIAG_DF_LOC <- file.path(DATA_DIR, 'diagnostic.Rdata')
RESULTS_DIR <- file.path('results')
CLEAN_DIR = TRUE # SET TO TRUE TO DELETE the DATA_DIR

TPRED = 7 # Number of days forward to predict (from last data day, not from Sys.Date())
#SPL_K = 7 # Number of spline basis functions to use. basis has approx K shifts in sign
SPL_K = floor(as.numeric(as.Date(DATE)+TPRED-as.Date(DATE_0))/7)-1 # My logic here is we want fewer than one peak per week.

# DATE = '2020-04-04' # Date of run
# NYT_FILE <- '../data/2020-04-03-covid19-nyt.csv'
# ACS_FILE <- '../data/us-acs.RData'
# DATE_0 <- '2020-03-01' # First date to use
# # DATE_0 doesn't assume that all series go back to DATE_0,
# # but it will filter away anything before.
# 
# ROOT <- '/home/nnagle/Dropbox/students/Piburn/covid-model'
# SAMPLES_ROOT <- file.path(ROOT, 'tmp') # All samples will be stored under {SAMPLES_ROOT}/{DATE}/{STATE}
# DATA_DIR <- file.path(ROOT, 'tmp',DATE)
# DIAG_DF_LOC <- file.path(DATA_DIR, 'diagnostic.Rdata')
# 
# RESULTS_DIR <- file.path(ROOT, 'results')
# 
# CLEAN_DIR = FALSE # SET TO TRUE TO DELETE the DATA_DIR
# 
# NITER = 3000
# NTHIN =3
# NCHAINS = 4
# NNODES = 30 # Number of nodes in cluster

# Probably don't mess with WARMUP
# WARMUP = 500
# WARMUP = (NITER/NTHIN)/2




