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
}

data {
  int                    N;      //number of observations
  int                    I;      //number of counties
  int                    T;      //number of time steps
  int                    TPRED;
  int                    Ki;      // number of county-level covariates
  int                    Kt;      // number of county-level covariates
  int                    Kit;      // number of unit-level covariates
  int<lower=1,upper=T-2> spl_K;  // rank of spline
  int                    J1;     // number of groups in level 1 (csa?)
  int                    J2;     // number of groups in level 2 (meto/nonmetro?)
  int<lower=0>           Y[N] ;  //Cases [i t Y]
  int<lower=1, upper=I>  Y_i[N]; //county of obs n
  int<lower=1, upper=T>  Y_t[N]; //time step of obs n
  int<lower=1,upper=I*T> Y_lin_idx[N];
  vector<lower=0>[I]     Pop;    // Population
  int<lower=0,upper=J1>  group1[I]; // level 1 groups (use zero to indicate no group)
  int<lower=0,upper=I>   N1; // number of counties with a level 1
  int<lower=0,upper=I>   group1_bin[N1];
  int<lower=0,upper=I>   group1_bin_id[N1];
  int<lower=0,upper=J2>  group2[I]; // level 2 groups (use zero to indicate no group)
  matrix[I,Ki]           Xi;      // county-level covariates
  matrix[T+TPRED,Kt]     Xt;      // county-level covariates
  matrix[N*(T+TPRED),Kit]          Xit;      // county-level covariates
  matrix[T,spl_K]        Z_spl;  // Design matrix of spline basis (rotated so coefs are N(0,1))
  matrix[T,TPRED]        krig_wt;
  matrix[TPRED,TPRED]    krig_var_chol;
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
  vector[Ki]       theta_i;     // coefficients on county-level covariates.
  vector[Kt]       theta_t;     // coefficients on date-level covariates.
  vector[Kit]       theta_it;     // coefficients on unite-level covariates.
  real            a;         // overall intercept (also baseline of group2);
  vector[I]       a0;        // county-level intercepts
  vector[J1]      a1;        // group-1-level intercepts
  vector[J2-1]    raw_a2;        // group-2-level intercepts
  real<lower=0>   tau_a0;    // scale of county-level intercepts
  real<lower=0>   tau_a1;    // scale of group1-level intercepts
  real<lower=0>   tau_a2;    // scale of group2-level intercepts
  real<lower=0>   raw_tau_splb0;  // unscaled: sd on spline coefs
  real<lower=0>   raw_tau_splb1;  // unscaled: sd on spline coefs
  real<lower=0>   raw_tau_splb2;  // unscaled: sd on spline coefs
  real<lower=0>   phi;   // inv_phi is overdispersin
  vector[I]       log_phi_i;//
  matrix[I,spl_K]  spl_b0;     // county spline coefficients
  matrix[J1,spl_K] spl_b1;     // level 1 spline coefficients
  matrix[J2,spl_K] spl_b2;     // level 2 spline coefficients
}

transformed parameters{
  vector[J2] a2;
  real tau_splb0; // scale on Spline coefs
  real tau_splb1; // scale on Spline coefs
  real tau_splb2; // scale on Spline coefs
  vector[I]  phi_i;
  
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
  vector[N] log_mu;
  matrix[I,T] log_lambda; // log cases per person. Spline only. No day or county effects
  matrix[I,T] spl0;
  matrix[J1,T] spl1;
  matrix[J2,T] spl2;
  
  spl0 = rep_matrix(tau_a0*a0,T) + (tau_splb0*spl_b0) * Z_spl';
  log_lambda = a + spl0;
  
  if(J1>0) {
    spl1 = rep_matrix(tau_a1*a1,T) + (tau_splb1*spl_b1) * Z_spl';
    log_lambda[group1_bin,] += spl1[group1_bin_id,];
  }
  
  if(J2>0) {
    spl2 = rep_matrix(tau_a2*a2,T) + (tau_splb2*spl_b2) * Z_spl';
    log_lambda += spl2[group2,];
  }
  
  if(Ki>0){
    log_lambda += rep_matrix(Xi * theta_i, T);
  }
  
  if(Kt>0){
    log_lambda += rep_matrix(to_row_vector(Xt[1:T,]*theta_t), I);
  }
  
  log_mu = log_pop[Y_i] + to_vector(log_lambda)[Y_lin_idx];
  if(Kit>0){
    log_mu += Xit*theta_it;
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
  theta_i ~ normal(0, .5); // the range of +/- 1 ~ 1 order of magnitude
  theta_t ~ normal(0, .5); // the range of +/- 1 ~ 1 order of magnitude
  theta_it ~ normal(0, .5); // the range of +/- 1 ~ 1 order of magnitude
  
  /// Likelihood
  target += neg_binomial_2_log_lpmf(Y | log_mu, phi_i[Y_i] .* exp(log_mu));
  
  //target += nb1_zero_lpmf(Y | exp(log_mu), sum(exp(log_mu0)), phi);
}

 generated quantities{
  int Y_sim[sample_flag ? I : 0,
            sample_flag ? T+TPRED : 0];
  int Ypred[N];
  matrix[I,T+TPRED] log_lambda; // log cases per person. Spline only. No day effects
  matrix[I,T+TPRED] spl0;
  matrix[J1,T+TPRED] spl1;
  matrix[J2,T+TPRED] spl2;
  real mean_theta_t;
  
  vector[ lppd_flag ? N  : 0 ]  lppd;
  vector[ lppd_flag ? N  : 0 ]  log_mu_i;
  vector[ post_pred ? N  : 0 ]  Y_i_sim;

  // Output three splines: spline2, spline1, and log_lamdba
  // splines 1 and 2 are missing county-level covariates and will not be on per/person scale 
  
  {
    matrix[I, T] data_spl;
    matrix[I, TPRED] pred_spl;
    data_spl = rep_matrix(tau_a0*a0,T) + (tau_splb0*spl_b0) * Z_spl';
    pred_spl = data_spl * krig_wt;
    for(j in 1:I){
      pred_spl[j,] = to_row_vector(multi_normal_cholesky_rng(to_vector(pred_spl[j,]), tau_splb0*krig_var_chol));
    }
    spl0 = append_col(data_spl, pred_spl);
  }
  log_lambda = a + spl0;
  
  if(J1>0){
    matrix[J1, T] data_spl;
    matrix[J1, TPRED] pred_spl;
    data_spl = rep_matrix(tau_a1*a1,T) + (tau_splb1*spl_b1) * Z_spl';
    pred_spl = data_spl * krig_wt;
    for(j in 1:J1){
      pred_spl[j,] = to_row_vector(multi_normal_cholesky_rng(to_vector(pred_spl[j,]), tau_splb1*krig_var_chol));
    }
    spl1 = append_col(data_spl, pred_spl);
    log_lambda[group1_bin,] += spl1[group1_bin_id,];
  }
  
  if(J2>0){
    matrix[J2, T] data_spl;
    matrix[J2, TPRED] pred_spl;
    data_spl = rep_matrix(tau_a2*a2,T) + (tau_splb2*spl_b2) * Z_spl';
    pred_spl = data_spl * krig_wt;
    for(j in 1:J2){
      pred_spl[j,] = to_row_vector(multi_normal_cholesky_rng(to_vector(pred_spl[j,]), tau_splb2*krig_var_chol));
    }
    spl2 = append_col(data_spl, pred_spl);
    log_lambda += spl2[group2,];
  }
  
  if(Ki>0){
    log_lambda += rep_matrix(Xi * theta_i, T+TPRED);
  }
  
  // Xt is added after sample_flag part
  // Xit not included yet
  
  if(sample_flag){
    for(i in 1:I){
      for(t in 1:(T+TPRED)){
        real est;
        real pe;
        int sim;
        est = min([log_pop[i] + log_lambda[i,t], 10.]);
        sim = neg_binomial_2_log_rng(est, phi_i[i] * exp(est));
        Y_sim[i,t] = is_nan(sim) ? -1 : sim; 
      }
    }
  }

//. The following is useful for posterior checks in shinystan
  if(lppd_flag || post_pred){
     log_mu_i = log_pop[Y_i] + to_vector(log_lambda)[Y_lin_idx] + Xt[Y_i] * theta_t;
     for(n in 1:N){
       real log_mu;
       real est;
       est = min([log_mu_i[n], 10.]);
       if(lppd_flag) lppd[n] = neg_binomial_2_log_lpmf(Y[n] | est, phi_i[Y_i[n]] * exp(est));
       if(post_pred){
         int sim;
         sim = neg_binomial_2_log_rng(est,phi_i[Y_i[n]] * exp(est));
         Y_i_sim[n] = is_nan(sim) ? -1 : sim; 
       } 
     }
  }
  
  mean_theta_t = (Kt>0) ? mean(exp(Xt*theta_t)) : 1;
  log_lambda += log(mean_theta_t);


}
