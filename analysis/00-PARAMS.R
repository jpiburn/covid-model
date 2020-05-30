DATE <- '2020-05-20' # Date of run
TN_ONLY <- TRUE
CLEAN_DIR = TRUE # SET TO TRUE TO DELETE the DATA_DIR
NYT_FILE <- NULL
ACS_FILE <- NULL
DATE_0 <- '2020-03-01' # First date to use
##ROOT = ""
SAMPLES_ROOT <- ifelse(Sys.info()[["nodename"]] == 'quetelet',
                       '/data/covid/tmp',
                       'tmp') # All samples will be stored under {SAMPLES_ROOT}/{DATE}/{STATE}
if(TN_ONLY) SAMPLES_ROOT <- file.path(SAMPLES_ROOT, 'TN')
NITER = 3000
NTHIN = 2
NCHAINS = 4
WARMUP = (NITER/NTHIN)/2
WARMUP = 500 # WARMUP is number of retained (after thinning) samples
NNODES = ifelse(Sys.info()[["nodename"]] == 'quetelet',
                30,
                80) # Number of nodes in cluster
if(TN_ONLY) NNODES = NCHAINS
DATA_DIR <- file.path(SAMPLES_ROOT,DATE)
DIAG_DF_LOC <- file.path(DATA_DIR, 'diagnostic.Rdata')
RESULTS_DIR <- file.path('results')

ZERO_PAD = 7 # Number of days of zeros to add before first obs for each county
TPRED = 7 # Number of days forward to predict (from last data day, not from Sys.Date())
#SPL_K = 7 # Number of spline basis functions to use. basis has approx K shifts in slope (approx K/2 peaks and K/2 valleys)
SPL_K = floor(as.numeric(as.Date(DATE)+TPRED-as.Date(DATE_0))/7)-1 # My logic here is we want fewer than one peak per week.
#SPL_K = floor(as.numeric(as.Date(DATE)+TPRED-as.Date(DATE_0))/4)-1 # Trying for a little more flexibility than above.






