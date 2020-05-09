Covid Model 1
================
Nicholas Nagle
3/25/2020

# Notation:

Data:

\(Y_{it}\) number of cases in county \(i\) on data \(t\).

\(P_i\) population of county \(i\).

\(X_i\): k-vector of covariates for county \(i\).

Latent rates:

\(\operatorname{Exp}[Y_{it}] = P_i \lambda_{it}\)

\(\lambda_{it}\): rate (per person) in county \(i\) on date \(t\)

## Basic idea

The count of cases is a poisson rv with some space-time varying rate.

For an exponentially growing disease, we model the rate as log-linear in
time.

“Flattening the curve” can be detected by a declining slope in
log-linear space.

So we’ll use:

\(\log(\lambda_{i,t}) = a_{i,t} + b_{i,t} t + \theta X_i\)

\(a_{it} = a_{i,t-1} + \epsilon^{[a]}_{i,t}\)

\(b_{it} = b_{i,t-1} + \epsilon^{[b]}_{i,t}\)

We are especially interested in places/times where
\(\epsilon^{[b]}_{i,t}\) is negative, as that represents “flattening”
the curve.

The following model has two levels for \(a\) and \(b\) and adds
correlation between \(\epsilon^{[a]}_{i,t}\) and
\(\epsilon^{[b]}_{i,t}\).

``` r
library(tidyverse)
library(rstan)
library(tidybayes)
```

``` stan
data{
  int N; //number of counties
  int T; //number of time steps
  int K; // number of county-level covariates
  int J1; // number of groups in level 1 (metros?)
  int Y[N,T]; //Cases
  vector[N] Pop; // Population
  int<lower=1,upper=J1> group1[N]; // level 1 groups
  matrix[N,K] X; // county-level covariates
}
transformed data{
  vector[N] log_pop;
  matrix[T,2] X_time;
    for(t in 1:T){
      X_time[t,1] = 1;
      X_time[t,2] = t-1;
    }
  log_pop = log(Pop);
}
parameters{
  vector[2] eps[N,T-1]; // level0 time series difference
  vector[2] eps_1[J1,T-1]; // level1 time series difference
  vector[2] ab0[N]; // level 0 time 0 cofficient
  vector[2] ab0_1[N]; /// level 1 time 0 coefficient
  vector[K] theta; // coefficients on county-level covariates.
  cholesky_factor_corr[2] L_eps; // correlation of eps
  cholesky_factor_corr[2] L_eps_1;
  vector<lower=0>[2] sigma_eps;
  vector<lower=0>[2] sigma_eps_1;
}
transformed parameters{
  vector[2] ab[N,T];
  vector[N] Xtheta;
  matrix[N,T] log_lambda;
  Xtheta = X*theta;
  for(i in 1:N){
    ab[i,1] = ab0_1[group1[i]];
    for(t in 2:T){
      ab[i,t] = ab[i,t-1] + sigma_eps_1 .* (L_eps_1 * eps_1[group1[i],t]) + 
                            sigma_eps   .* (L_eps   * eps[i,t]);
    }
  }
  for(i in 1:N){
    for(t in 1:T){
      log_lambda[i,t] = ab[i,t,1] + ab[i,t,2]*t + Xtheta[i];
    }
  }
}
model{
  theta ~ normal(0,10);
  L_eps ~ lkj_corr_cholesky(2.0);
  L_eps_1 ~ lkj_corr_cholesky(2.0);
  sigma_eps ~ cauchy(0., 1.);
  sigma_eps_1 ~ cauchy(0., 1.);
  for(i in 1:N){
    ab0[i] ~ normal(0,10);
    for(t in 1:(T-1)) eps[i,t] ~ normal(0,1);
  }
  for(j in 1:J1){
    ab0_1[j] ~ normal(0,10);
    for(t in 1:(T-1)) eps_1[j,t] ~ normal(0,1);
  }
  for(i in 1:N) {
    for(t in 2:T){
      Y[i,t] ~ poisson_log(log_pop[i] + log_lambda[i,t]);
    }
  }
}
```
