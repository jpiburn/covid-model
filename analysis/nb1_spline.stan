functions{
  real nb1_log_lpmf(int Y, real log_mu, real phi){
    // create NB1 wih variance = mu+ (1/phi)* mu
    real logp;
    real scale = phi*exp(log_mu);
    //logp = lchoose(Y + scale - 1, Y);
    logp = lgamma(Y + scale) - lgamma(scale);
    logp += -Y*log1p(phi) - scale*log1p(1/phi);
    return(logp);
  }
  real nb1_zero_lpmf(int[] y, real[] mu, real sum_mu0, real k){
    // create NB1 wih variance = mu+ (1/phi)* mu
    // add the 0 entries separately to simplify comput
    // vector of y (nonzero y's)
    // vector of mu (for nonzero y's)
    // sum of mu (for zero y's)
    int N = size(y);
    //vector[N] phi;
    real lp = 0;
    //real sum_y = sum(y);
    ///real sum_phi ;
    //real sum_phi0 = sum_mu0*k;
    //for(i in 1:N) phi[i] = mu[i]*k;
    //sum_phi = sum(phi);
    // Term owing to y>1
    for( i in 1:N) {
      real scale = mu[i]*k;
      lp += lgamma(y[i] + scale) - lgamma(scale) - lgamma(y[i]+1);
      //lp += lchoose(y[i] + phi[i] -1, y[i]);
      lp += -y[i]*log1p(k) - scale*(log1p(k) - log(k));
    }
    //lp += -sum_y*log1p(k) - sum_phi*(log1p(k) - log(k));
    // Term owing to y=0
    lp += - k* sum_mu0*log1p(1/k);
    return(lp);
  }
}

data {
  int                    N;      //number of observations
  int                    N0;     // number of 0 counts
  int                    I;      //number of counties
  int                    T;      //number of time steps
  int                    K;      // number of county-level covariates
  int<lower=1,upper=T-2> spl_K;  // rank of spline
  int                    J1;     // number of groups in level 1 (metros?)
  int                    J2;     // number of groups in level 2 (csas?)
  int<lower=0>           Y[N] ;  //Cases [i t Y]
  int<lower=1, upper=I>  Y_i[N]; //county of obs n
  int<lower=1, upper=T>  Y_t[N]; //time step of obs n
  int<lower=1, upper=I>  Y0_i[N0]; 
  int<lower=1, upper=T>  Y0_t[N0]; 
  vector<lower=0>[I]     Pop;    // Population
  int<lower=0,upper=J1>  group1[I]; // level 1 groups (use zero to indicate no group)
  int<lower=0,upper=J2>  group2[I]; // level 2 groups (use zero to indicate no group)
  matrix[I,K]            X;      // county-level covariates
  matrix[T,spl_K]        Z_spl;  // Design matrix of spline basis (rotated so coefs are N(0,1))
  matrix[T,7]            X_dow;  // Day of week design matrix
  real<lower=0>          taub0_scale;
  real<lower=0>          taub1_scale;
  real<lower=0>          taub2_scale;
  int                    sample_flag; //set to 0 if you don't want samples
  int<lower=0, upper=1>  lppd_flag; // set to 0 if you don't want log posterior predictive
  int<lower=0, upper=1>  post_pred; // set to 0 if you don't want posterior sim
}

transformed data{
  vector[I]     log_pop;
  
  log_pop = log(Pop);
}

parameters{
  vector[K]       theta;     // coefficients on county-level covariates.
  real            a;         // overall intercept;
  vector[I]       a0;        // county-level intercepts
  vector[J1]      a1;        // group-1-level intercepts
  vector[J2-1]    raw_a2;        // group-2-level intercepts
  real<lower=0>   tau_a0;    // scale of county-level intercepts
  real<lower=0>   tau_a1;    // scale of group1-level intercepts
  real<lower=0>   tau_a2;    // scale of group2-level intercepts
  real<lower=0>   raw_tau_splb0;  // unscaled: sd on spline coefs
  real<lower=0>   raw_tau_splb1;  // unscaled: sd on spline coefs
  real<lower=0>   raw_tau_splb2;  // unscaled: sd on spline coefs
  real<lower=0>   phi;   // inv_phi is on scale of glm
  vector[I]       log_phi_i;//
  matrix[I,spl_K]  spl_b0;     // county spline coefficients
  matrix[J1,spl_K] spl_b1;     // level 1 spline coefficients
  matrix[J2,spl_K] spl_b2;     // level 2 spline coefficients
  //real<lower=.01> a_gamma;   // inverse of Neg Bin scale (there is no way that phi>100)
}

transformed parameters{
  vector[J2] a2;
  real tau_splb0; // scale on Spline coefs
  real tau_splb1; // scale on Spline coefs
  real tau_splb2; // scale on Spline coefs
  real inv_phi;
  vector[I]  phi_i;
  
  //real phi;               // Neg Bin scale
  inv_phi = 1/phi;
  //phi = 1/a_gamma;        // var NB(2) = mu + mu^2/phi
                          // var NB(1) = mu + mu/phi
  a2[1] = 0;
  if(J2>1){
    for(j in 2:J2){
      a2[j] = raw_a2[j-1];
    }
  }
  for(i in 1:I){
    phi_i[i] = exp(.2*log_phi_i[i])*phi;
  }
  tau_splb0 = taub0_scale * raw_tau_splb0;
  tau_splb1 = taub1_scale * raw_tau_splb1;
  tau_splb2 = taub2_scale * raw_tau_splb2;

}


model {
  real log_mu[N];
  real log_mu0[N0];
  matrix[I,T] log_lambda; // log cases per person. Spline only. No day or county effects
  for(i in 1:I){
    log_lambda[i]  =  a + tau_a0 * a0[i] + 
                      tau_splb0 * to_row_vector(Z_spl * to_vector(spl_b0[i]));
    if(J1>0){ // Add csa effects if present
      if(group1[i] > 0){ 
        log_lambda[i] += tau_a1 * a1[group1[i]] + 
                         tau_splb1 * to_row_vector(Z_spl * to_vector(spl_b1[group1[i]]));
      }
    }
    if(J2>0){ // Add metro effects if present
      if(group2[i] > 0){
        log_lambda[i] += tau_a2*a2[group2[i]] + 
                         tau_splb2 * to_row_vector(Z_spl * to_vector(spl_b2[group2[i]]));
      }
    }  
  }
  
  a ~ normal(-11,10);
  raw_a2 ~ normal(0,1);  // exp(-11) ~ .1 cases in 10,000 persons we're going to put the overall intercept in the metro level
  a1 ~ normal(0,1);
  a0 ~ normal(0,1);
  
  to_vector(spl_b0) ~ normal(0,1);
  to_vector(spl_b1) ~ normal(0,1);
  to_vector(spl_b2) ~ normal(0,1);

  // Prior for intercepts of each level
  tau_a0 ~ normal(0.,1.5); //+/- 4 -> [.018-54] times.  Allows ~ 3 orders of magnitude between large and small
  tau_a1 ~ normal(0.,1.5);
  tau_a2 ~ normal(0.,1.5); 

  // prior for spline scale
  raw_tau_splb0 ~ normal(0,1);
  raw_tau_splb1 ~ normal(0,1);
  raw_tau_splb2 ~ normal(0,1);
  
  phi ~ normal(.25,.5); // DON'T FORGET: phi in inverse relative to glm()
  log_phi_i ~ normal(0,1); // N(0,.2): most of exp(log_phi_i) will be between .5 and 2 
  theta ~ normal(0, .5); // the range of +/- 1 ~ 1 order of magnitude
  
  /// Likelihood
  // likelihood term for non-zero counts
  for(n in 1:N) {
    log_mu[n] = log_pop[Y_i[n]] + log_lambda[Y_i[n],Y_t[n]] +
                              X[Y_i[n]] * theta;
    target += neg_binomial_2_log_lpmf(Y[n] | log_mu[n], phi_i[Y_i[n]] * exp(log_mu[n]));
    //target += -Y[n]*log1p(phi) - phi*mu*(log1p(phi) - log(phi));
  }
  // likelihood term for counts of zero
  for(n0 in 1:N0) {
    log_mu0[n0] = log_pop[Y0_i[n0]] + log_lambda[Y0_i[n0],Y0_t[n0]] +
              X[Y0_i[n0]] * theta;  
    target += neg_binomial_2_log_lpmf(0 | log_mu0[n0], phi_i[Y0_i[n0]] * exp(log_mu0[n0]));
    //target += neg_binomial_2_log_lpmf(0 | log_mu0, phi*exp(log_mu0));
    //target += -phi*exp(log_mu0)*log1p(1/phi);
  }
  //target += nb1_zero_lpmf(Y | exp(log_mu), sum(exp(log_mu0)), phi);
}

 generated quantities{
  int Y_sim[sample_flag ? I : 0,
            sample_flag ? T : 0];
  int Ypred[N];
  matrix[I,T] log_lambda; // log cases per person. Spline only. No day effects
  matrix[J1,T] spl1;
  matrix[J2,T] spl2;
  vector[ lppd_flag ? N  : 0 ]  lppd;
  vector[ lppd_flag ? N0 : 0 ]  lppd0;
  vector[ lppd_flag ? N  : 0 ]  log2_mu;
  vector[ lppd_flag ? N0 : 0 ]  log2_mu0;

  vector[ post_pred ? N  : 0 ]  Yi_sim;
  vector[ post_pred ? N0 : 0 ]  Yi0_sim;
  
  // Output three splines: spline2, spline1, and log_lamdba
  // splines 1 and 2 are missing county-level covariates and will not be on per/person scale 
  if(J1>0){
    for(j in 1:J1){
      spl1[j] = tau_a1 * a1[j] + 
                tau_splb1 * to_row_vector(Z_spl * to_vector(spl_b1[j]));
    }
  }
  if(J2>0){
    for(j in 1:J2){
      spl2[j] = tau_a2 * a2[j] + 
                tau_splb2 * to_row_vector(Z_spl * to_vector(spl_b2[j]));
    }
  }
  for(i in 1:I){
    log_lambda[i] =  a + tau_a0 * a0[i] + 
                     tau_splb0 * to_row_vector(Z_spl * to_vector(spl_b0[i]));
    if(J1>0){
      if(group1[i] > 0){
        log_lambda[i] += spl1[group1[i]];
      }
    } 
    if(J2>0){
      if(group2[i] > 0){
        log_lambda[i] += spl2[group2[i]];
      }
    } 
    log_lambda[i] += X[i]*theta;
    ///+ log(mean(exp(theta_dow)));
  }
                     
  if(sample_flag){
    for(i in 1:I){
      for(t in 1:T){
        real log_mu;
        real est;
        real pe;
        int sim;
        log_mu = log_pop[i] + log_lambda[i,t] ;
        est = min([log_mu, 10.]);
        sim = neg_binomial_2_log_rng(est, phi_i[i] * exp(est));
        Y_sim[i,t] = is_nan(sim) ? -1 : sim; 
      }
    }
  }

//. The following is useful for posterior checks in shinystan
  if(lppd_flag || post_pred){
     for(n in 1:N){
       real log_mu;
       real est;
       log2_mu[n] = log_pop[Y_i[n]] + log_lambda[Y_i[n],Y_t[n]] +
                    X[Y_i[n]] * theta;
       est = min([log2_mu[n], 10.]);
       if(lppd_flag) lppd[n] = neg_binomial_2_log_lpmf(Y[n] | est, phi_i[Y_i[n]] * exp(est));
       if(post_pred){
         int sim;
         sim = neg_binomial_2_log_rng(est,phi_i[Y_i[n]] * exp(est));
         Yi_sim[n] = is_nan(sim) ? -1 : sim; 
       } 
     }
     for(n in 1:N0){
       real log_mu;
       real est;
       log2_mu0[n] = log_pop[Y0_i[n]] + log_lambda[Y0_i[n],Y0_t[n]] +
                     X[Y0_i[n]] * theta;
       est = min([log2_mu0[n], 10.]);
       if(lppd_flag) lppd0[n] = neg_binomial_2_log_lpmf(0 | est, phi_i[Y0_i[n]] * exp(est));
       if(post_pred) Yi0_sim[n] = neg_binomial_2_log_rng(est, phi_i[Y0_i[n]] * exp(est));
     }
  }

}
