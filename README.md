
<!-- README.md is generated from README.Rmd. Please edit that file -->

# covid-model

<!-- badges: start -->

<!-- badges: end -->

An initial repo to test spatio-temporal modeling of covid-19.

## Example Data

There are two files under `data/` relating to Washington state

  - `washington-acs.RData`: data.frame including `sf` geometries.
    Contains example variables from ACS as well as county membership of
    MSAs
  - `washington-covid19.csv`: csv of daily county level covid19 cases
    and deaths.

Both of these files can be recreated with the scripts in `data-raw/`
