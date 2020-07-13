# This is used b/c of weird permission issues when running the sricpt through Rstudio vs R cmd
# I run this script as a background job, but if the write permissions issue were fixed
# then running analysis/01-run-full-model-pipeline.R directly would be what to do

system('sudo Rscript analysis/01-run-full-model-pipeline.R')
