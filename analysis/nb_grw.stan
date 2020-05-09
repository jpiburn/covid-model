data{
  int N; //number of observations
  int I; //number of counties
  int T; //number of time steps
  int K; // number of county-level covariates
  int J1; // number of groups in level 1 (metros?)
  int<lower=0> Y[N] ; //Cases [i t Y]
  int<lower=1, upper=I> Y_i[N]; //county of obs n
  int<lower=1, upper=T> Y_t[N]; //time step of obs n
  vector<lower=0>[I] Pop; // Population
  int<lower=0,upper=J1> group1[I]; // level 1 groups
  matrix[I,K] X; // county-level covariates
  real<lower=0> taub_scale;
}
transformed data{
  vector[I] log_pop;
  row_vector[T] t_vec;
  matrix[T-1,T-1] X_time; // piecewise spline matrix
  log_pop = log(Pop);
  for(r in 1:(T-1)){
    for(c in 1:(T-1)){
      X_time[r,c] = (r>=c) ? r-c + 1 : 0; // lower diagonal toeplitz: 1 2 3 ... 
    }
  }
  for(t in 1:T) t_vec[t] = t-1;
}
parameters{
  vector[K] theta; // coefficients on county-level covariates.
  real a; // global intercept
  real b; // global slope
  vector<upper=3>[I] a0; // county-level intercepts
  vector<upper=3>[J1] a1; // group-1-level intercepts
  real<lower=0,upper=3> tau0; // scale of county-level intercepts
  real<lower=0,upper=3> tau1; // scale of group1-level intercepts
  real<lower=0,upper=.3> taub; // scale on spline coefficients
  matrix[I,T-1] b0_raw; // spline coefficients
  real<lower=.01> a_gamma; // inverse of Neg Bin scale (there is no way that phi>100)
}
transformed parameters{
  vector[I] Xtheta;
  matrix[I,T] log_lambda; // log cases per person
  real phi; // Neg Bin scale
  Xtheta = X*theta;
  phi = 1/a_gamma;
  for(i in 1:I){
    log_lambda[i] = Xtheta[i] + a + tau0*a0[i] + b*t_vec;
    if(J1>0){ // Check for states with no metros.
      log_lambda[i] += tau1*a1[group1[i]];
      }
    log_lambda[i,2:] += taub * to_row_vector(X_time * to_vector(b0_raw[i]));
  }
}
model{
  a ~ normal(-5,10);
  a0 ~ normal(0.,1.);
  a1 ~ normal(0., 1.);
  b ~ normal(0,1);
  tau0 ~ cauchy(0.,.1);
  tau1 ~ cauchy(0.,1.);
  taub ~ cauchy(0, taub_scale);
  to_vector(b0_raw) ~ student_t(3,0.,1.);
  //a_gamma ~ gamma(.001,.001);
  a_gamma ~ cauchy(0,1.);
  for(n in 1:N) {
    Y[n] ~ neg_binomial_2_log(log_pop[Y_i[n]] + log_lambda[Y_i[n],Y_t[n]],phi);
  }
}
generated quantities{
 matrix[I,T] log_mu;
 matrix[I,T-1] cum_slope;
 int Y_sim[I,T];
 for(i in 1:I){
   cum_slope[i] = taub * cumulative_sum(b0_raw[i]);
   for(t in 1:T){
     real est;
     log_mu[i,t] = log_pop[i]+log_lambda[i,t];
     est = min([log_mu[i,t], 15.]);
     Y_sim[i,t] = neg_binomial_2_log_rng(est, phi);
   }
 }
}
