functions{
  real nb1_log_lpmf(int Y, real log_mu, real phi){
    real logp;
    real scale = phi*exp(log_mu);
    logp = lchoose(Y + scale - 1, Y);
    logp += -Y*log1p(phi) - scale*log1p(1/phi);
    return(logp);
  }
}

data {
  int                    N;      //number of observations
  int                    I;      //number of counties
  int                    T;      //number of time steps
  int                    K;      // number of county-level covariates
  int<lower=1,upper=T-2> spl_K;  // rank of spline
  int                    J1;     // number of groups in level 1 (metros?)
  int<lower=0>           Y[N] ;  //Cases [i t Y]
  int<lower=1, upper=I>  Y_i[N]; //county of obs n
  int<lower=1, upper=T>  Y_t[N]; //time step of obs n
  vector<lower=0>[I]     Pop;    // Population
  int<lower=0,upper=J1>  group1[I]; // level 1 groups
  matrix[I,K]            X;      // county-level covariates
  matrix[T,spl_K]        Z_spl;  // Design matrix of spline basis (rotated so coefs are N(0,1))
  matrix[T,7]            X_dow;  // Day of week design matrix
  real<lower=0>          taub_scale;
}

transformed data{
  vector[I]     log_pop;
  row_vector[T] t_vec; // Just a vector from 0:T
  //real phi = 0.05;
  
  log_pop = log(Pop);
  for(t in 1:T) t_vec[t] = t-1;
}

parameters{
  vector[K]       theta;     // coefficients on county-level covariates.
  vector[7]       theta_dow; // coefficients on day of week
  real            a;         // global intercept
  real            b;         // global slope
  vector[I]       a0;        // county-level intercepts
  vector[J1]      a1;        // group-1-level intercepts
  vector[I]       b0;        // county-level intercepts
  real<lower=0>   tau_a0;    // scale of county-level intercepts
  real<lower=0>   tau_a1;    // scale of group1-level intercepts
  real<lower=0>   raw_tau_b0;    // scale of county-level slopes
  real<lower=0>   raw_tau_splb;  // scale on spline coefficients
  real<lower=0>   phi;   // inv_phi is on scale of glm
  matrix[I,spl_K] spl_b;     // spline coefficients
  //real<lower=.01> a_gamma;   // inverse of Neg Bin scale (there is no way that phi>100)
}

transformed parameters{
  real tau_splb;
  real tau_b0;
  real inv_phi;
  matrix[I,T] log_lambda; // log cases per person. Spline only. No day or county effects
  //real phi;               // Neg Bin scale
  inv_phi = 1/phi;
  //phi = 1/a_gamma;        // var NB(2) = mu + mu^2/phi
                          // var NB(1) = mu + mu/phi
  tau_splb = taub_scale * raw_tau_splb;
  tau_b0 = .04 * raw_tau_b0;
  for(i in 1:I){
    log_lambda[i] =  a + tau_a0*a0[i] + 
                    (b + tau_b0*b0[i]) * t_vec;
    log_lambda[i,] += tau_splb * to_row_vector(Z_spl * to_vector(spl_b[i]));
    if(J1>0){ // Add metro effect if metros are present
      log_lambda[i] += tau_a1*a1[group1[i]];
      }
  }
}


model {
  a ~ normal(-11,10);  // exp(-11) ~ .1 cases in 10,000 persons
  b ~ normal(.2, .1); // 3 orders of magnitude over 30 days = log(1000)/30 ~ .23
  
  a0 ~ normal(0,1);
  a1 ~ normal(0,1);
  b0 ~ normal(0,1);
  to_vector(spl_b) ~ normal(0,1);

  tau_a0 ~ normal(0.,1.5); //+/- 4 -> [.018-54] times.  Allows ~ 3 orders of magnitude between large and small
  tau_a1 ~ normal(0.,1.5);
  raw_tau_b0 ~ normal(0,1);
  phi ~ normal(.25,.5); // DON'T FORGET: phi in inverse relative to glm()
  //tau_b0 ~ normal(0, .08); // +/- .16 allows ~ +/-2 order of magnitude change over 30 days
  
  //tau_splb ~ normal(0, taub_scale);
  raw_tau_splb ~ normal(0, 1);
  //a_gamma ~ cauchy(0,1.); /// fat tails allow small scale
  
  theta ~ normal(0, .5); // the range of +/- 1 ~ 1 order of magnitude
  theta_dow ~ normal(0, .5); // the range of +/- 1 ~ 1 order of magnitude
  
  /// Likelihood
  for(n in 1:N) {
    real log_mu;
    log_mu = log_pop[Y_i[n]] + log_lambda[Y_i[n],Y_t[n]] +
                              X[Y_i[n]] * theta + 
                              X_dow[Y_t[n]] * theta_dow;
    //target += lchoose(Y[n] + phi*exp(log_mu) - 1, Y[n]);
    //target += -Y[n]*log1p(phi) - exp(log_mu)*phi*log1p(1/phi);
    Y[n] ~ neg_binomial_2_log(log_mu, phi * exp(log_mu) );
  }
}

generated quantities{
 int Y_sim[I,T];
 int Ypred[N];
 for(i in 1:I){
   for(t in 1:T){
     real log_mu;
     real est;
     int sim;
     log_mu = log_pop[i]+log_lambda[i,t] + 
              X[i]*theta +
              X_dow[t]*theta_dow;
     est = min([log_mu, 15.]);
     sim = neg_binomial_2_log_rng(est, phi * exp(log_mu));
     if (is_nan(sim)) 
       Y_sim[i,t] = -1;
     else 
       Y_sim[i,t] = sim;
   }
 }
 for(n in 1:N){
   real log_mu;
   int sim;
   real est;
   log_mu = log_pop[Y_i[n]] + log_lambda[Y_i[n],Y_t[n]] +
                              X[Y_i[n]] * theta + 
                              X_dow[Y_t[n]] * theta_dow;
   est = min([log_mu, 15.]);
   sim = neg_binomial_2_log_rng(est, phi * exp(est) );     
   if (is_nan(sim)) 
      Ypred[n] = -1;
   else 
      Ypred[n] = sim;
 }
}