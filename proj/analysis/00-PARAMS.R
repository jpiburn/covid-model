DATE <- '2020-06-24' # Date of run
TN_ONLY <- FALSE
NYT_FILE <- NULL
ACS_FILE <- NULL
DATE_0 <- '2020-03-01' # First date to use
DATE_N <- DATE # Last date to use 


SAMPLES_ROOT <- '/covidmodeldata'
NWORKERS <- 0
NNODES <- 32
RESULTS_DIR <- file.path('/covidmodeldata/results') 
CLEAN_DIR <- TRUE
SSH_KEY <- "/home/cades/covid-model/.ssh/id_rsa"
cades_workers <- 
  tibble::tribble(
    ~instance_name,    ~ip_address,
    "covid-worker-20", "172.22.3.223",
    "covid-worker-19", "172.22.3.218",
    "covid-worker-18", "172.22.3.222",
    "covid-worker-17", "172.22.3.219",
    "covid-worker-16", "172.22.3.212",
    "covid-worker-15", "172.22.3.217",
    "covid-worker-14", "172.22.3.220",
    "covid-worker-13", "172.22.3.211",
    "covid-worker-12", "172.22.3.209",
    "covid-worker-11", "172.22.3.221",
    "covid-worker-10", "172.22.3.216",
    "covid-worker-9", "172.22.3.214",
    "covid-worker-8", "172.22.3.208",
    "covid-worker-7", "172.22.3.210",
    "covid-worker-6", "172.22.3.215",
    "covid-worker-5", "172.22.3.207",
    "covid-worker-4", "172.22.3.206",
    "covid-worker-3", "172.22.3.205",
    "covid-worker-2", "172.22.3.204",
    "covid-worker-1", "172.22.3.203"
  )


NITER = 1500 # Number of  iterations (before thinning)
NTHIN = 2
NCHAINS = 8
WARMUP = (NITER/NTHIN)/2
WARMUP = 500 # WARMUP is number of retained (before thinning) samples

DATA_DIR <- file.path(SAMPLES_ROOT,DATE)
DIAG_DF_LOC <- file.path(DATA_DIR, 'diagnostic.Rdata')

ZERO_PAD = 7 # Number of days of zeros to add before first obs for each county
TPRED = 7 # Number of days forward to predict (from last data day, not from Sys.Date())
SPL_K = floor(as.numeric(as.Date(DATE)-as.Date(DATE_0))/14) # My logic here is we want fewer than one peak per week.




