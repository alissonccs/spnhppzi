#include /chunks/mylib.stan
data{
  int <lower=0> N;
  int <lower=0> p;
  int <lower=0>n;
  int <lower=0> baseline;
  int m;
  int id[N];
  int n_ind [n];
  int begin_ind [n];
  int end_ind[n];
  vector [N] event;
  vector [n] max_stop;
  vector [N] time;
  matrix [p == 0 ? 0 : N, p] X;
  matrix [p == 0 ? 0 : n, p] Xy;
  vector [n] IndRec2;
  int<lower=0> approach;
  int <lower=0> ZI;
  real shp_alpha1;
  real scl_alpha1;
  real shp_alpha2;
  real scl_alpha2;
  real mu_omega;
  real shp_sigma2_z;
  real scl_sigma2_z;
  real mu_beta;
  real <lower=0> sigma_beta;
  }

parameters{
  vector  <lower=0> [m] alpha;
  vector [p] beta;
  real <lower=0,upper=1> pii [ZI == 0 ? 0 : 1];
  vector[n] omega_1;
  real <lower=0> sigma2_z;
          }

 model {
  vector [n] Lambda0 ;
  vector [N] log_lambda0 ;
  vector [N] log_lambda0_event ;
  vector[p == 0 ? 0 : N] eta;
  vector [N] eta_event ;
  vector[p == 0 ? 0 : n] exp_etay;
  vector [n] sum_log_lambda0 = rep_vector(0, n);
  vector [n] sum_eta = rep_vector(0, n);
  int a = 0;
  int c = 0;

  if(p>0){
  for (i in 1:N){
     eta[i] = X[i,]*beta+omega_1[id[i]];
     eta_event[i] = event[i]*eta[i];
  }
  for (j in 1:n){
     exp_etay[j] = exp(Xy[j,]*beta+omega_1[j]);
        }
}


if(baseline==1){
   Lambda0 = Lambda_plp(max_stop, alpha,n);
    log_lambda0 = log_lambda_plp(time, N, alpha);
    log_lambda0_event = event .*log_lambda0;
}


 if(p == 0){
    for (i in 1:n) {

       if(IndRec2[i] == 0)
         target += log_sum_exp(bernoulli_lpmf(1| pii),
                  Lambda0[i] +
                   bernoulli_lpmf(0 | pii));
       else
         target += bernoulli_lpmf(0 | pii)+
                   Lambda0[i] +
                   sum(log_lambda0_event[begin_ind[i]:end_ind[i]]);
                  }
          }
    else{
      for (i in 1:n) {
       if(IndRec2[i] == 0)
         target += log_sum_exp(bernoulli_lpmf(1| pii),
                   Lambda0[i]*exp_etay[i] +
                   bernoulli_lpmf(0 | pii));
       else
          target +=bernoulli_lpmf(0 | pii)+
                  Lambda0[i]*exp_etay[i] +
                  sum(log_lambda0_event[begin_ind[i]:end_ind[i]])+
                  sum(eta_event[begin_ind[i]:end_ind[i]]);
                     }
        }



  if (approach==1){
            alpha[1] ~ gamma(shp_alpha1,scl_alpha1);
            alpha[2] ~ gamma(shp_alpha2,scl_alpha2);
            beta ~ normal(mu_beta,sigma_beta);
            sigma2_z ~ gamma(shp_sigma2_z,scl_sigma2_z);
            omega_1 ~ normal(mu_omega,sigma2_z);
            }
                                           }

generated quantities{
  vector[n] log_lik;
  {
    vector [n] Lambda0 ;
  vector [N] log_lambda0 ;
  vector [N] log_lambda0_event ;
  vector[p == 0 ? 0 : N] eta;
  vector [N] eta_event ;
  vector[p == 0 ? 0 : n] exp_etay;
  vector [n] sum_log_lambda0 = rep_vector(0, n);
  vector [n] sum_eta = rep_vector(0, n);
  int a = 0;
  int c = 0;


  if(p>0){
  for (i in 1:N){
     eta[i] = X[i,]*beta+omega_1[id[i]];
     eta_event[i] = event[i]*eta[i];
  }
  for (j in 1:n){
     exp_etay[j] = exp(Xy[j,]*beta+omega_1[j]);
        }
}

if(baseline==1){
   Lambda0 = Lambda_plp(max_stop, alpha,n);
    log_lambda0 = log_lambda_plp(time, N, alpha);
    log_lambda0_event = event .*log_lambda0;
}


 if(p == 0){
    for (i in 1:n) {

       if(IndRec2[i] == 0){
         log_lik[i]= log_sum_exp(bernoulli_lpmf(1| pii),
                  Lambda0[i] +
                   bernoulli_lpmf(0 | pii));}
       else{
         log_lik[i]= bernoulli_lpmf(0 | pii)+
                   Lambda0[i] +
                   sum(log_lambda0_event[begin_ind[i]:end_ind[i]]);
                  }
                  }
          }
    else{
      for (i in 1:n) {
       if(IndRec2[i] == 0){
         log_lik[i]= log_sum_exp(bernoulli_lpmf(1| pii),
                   Lambda0[i]*exp_etay[i] +
                   bernoulli_lpmf(0 | pii));
                   }
       else{
          log_lik[i]=bernoulli_lpmf(0 | pii)+
                  Lambda0[i]*exp_etay[i] +
                  sum(log_lambda0_event[begin_ind[i]:end_ind[i]])+
                  sum(eta_event[begin_ind[i]:end_ind[i]]);}
                     }
        }
  }
}
