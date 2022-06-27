#include /chunks/mylib.stan
data{
  int <lower=0> N;
  int <lower=0> p;
  int <lower=0>n;
  int <lower=0> tp_hf;
  int m;
  int id[N];
  int n_ind [n];
  int begin_int [n];
  int end_int[n];
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
  // real <lower=0> sigma_omega;
  real shp_sigma_omega;
  real scl_sigma_omega;
  real mu_beta;
  real <lower=0> sigma_beta;
  int <lower=0> tp_prior;
    }

parameters{
  vector  <lower=0> [m] alpha;
  vector [p] beta;
  real <lower=0,upper=1> pii [ZI == 0 ? 0 : 1];
  vector [n] omega;
  real <lower=0> sigma_omega;
          }


transformed parameters {
  vector[n] log_lik1;
  {
  vector [n] Lambda0 ;
  vector [N] log_lambda0 ;
  vector [N] log_lambda0_ ;
  vector[p == 0 ? 0 : N] eta;
  vector [N] eta_ ;
  vector[p == 0 ? 0 : n] exp_etay;
  vector [n] sum_log_lambda0 = rep_vector(0, n);
  vector [n] sum_eta = rep_vector(0, n);
  int a = 0;
  int c = 0;

if(p>0){
  for (i in 1:N){
     eta[i] = X[i,]*beta+omega[id[i]];
     eta_[i] = event[i]*eta[i];
  }
     exp_etay = exp(Xy*beta+omega);
        }

  Lambda0 = Lambda_plp2(max_stop, alpha,n);
  log_lambda0 = lambda_plp2(time, N, alpha);
  log_lambda0_ = event .*log_lambda0;


 // CALCULA VEROSSIMILHANÇA ACUMULADA POR INDIVÍDUO
 // for (b in 1:n) {
 //    sum_log_lambda0[b] = 0;
 //    sum_eta[b] = 0;
 //  for (i in 1:n_ind[b]) {
 //    if (a <= N) {
 //    a = a + 1;
 //    sum_log_lambda0[b] += log_lambda0_[a];
 //    if(p>0){
 //    sum_eta[b] += eta_[a];
 //    }
 //                 }
 //                         }
 //  }

 // for ( b in 1:n) {
 //       vector [n_ind[b]] sub_lambda = log_lambda0_[sum(n_ind[:b-1])+1:sum(n_ind[:b])];
 //       vector [n_ind[b]] sub_eta = eta_[sum(n_ind[:b-1])+1:sum(n_ind[:b])];
 //       sum_log_lambda0[b] = 0;
 //       sum_eta[b] = 0;
 //       sum_log_lambda0[b]=sum(log_lambda0_[sum(n_ind[:b-1])+1:sum(n_ind[:b])]);
 //       if(p>0){
 //       sum_eta[b]=sum(sub_eta);
 //               }
 //                 }


  for ( b in 1:n) {
        sum_log_lambda0[b]=sum(log_lambda0_[begin_int[b]:end_int[b]]);
        if(p>0){
       sum_eta[b]=sum(eta_[begin_int[b]:end_int[b]]);
               }
 }


 if(p == 0){
    for (i in 1:n) {

       if(IndRec2[i] == 0)
         log_lik1[i]= log_sum_exp(bernoulli_lpmf(1| pii),
                  Lambda0[i] +
                   bernoulli_lpmf(0 | pii));
       else
         log_lik1[i]= bernoulli_lpmf(0 | pii)+
                   Lambda0[i] +
                   sum_log_lambda0[i];
                  }
          }
    else{
      for (i in 1:n) {
       if(IndRec2[i] == 0)
         log_lik1[i]= log_sum_exp(bernoulli_lpmf(1| pii),
                   Lambda0[i]*exp_etay[i] +
                   bernoulli_lpmf(0 | pii));
       else
         log_lik1[i]=bernoulli_lpmf(0 | pii)+
                  Lambda0[i]*exp_etay[i] +
                  sum_log_lambda0[i]+sum_eta[i];
                     }
        }

 // if(approach==1 && tp_prior==1){
 //            alpha[1] ~ gamma(shp_alpha1,scl_alpha1);
 //            alpha[2] ~ gamma(shp_alpha2,scl_alpha2);
 //            beta ~ normal(mu_beta,sigma_beta);
 //            // sigma_omega ~ gamma(shp_sigma_omega,scl_sigma_omega);
 //            omega[id] ~ normal(mu_omega,sigma_omega);
 //                               }
  }
}

model{
     target +=log_lik1;

if(approach==1 && tp_prior==1){
            alpha[1] ~ gamma(shp_alpha1,scl_alpha1);
            alpha[2] ~ gamma(shp_alpha2,scl_alpha2);
            beta ~ normal(mu_beta,sigma_beta);
            sigma_omega ~ gamma(shp_sigma_omega,scl_sigma_omega);
            omega~ normal(-(sigma_omega)^2/2,sigma_omega);
            // omega ~ normal(mu_omega,sigma_omega);
                               }
  }
