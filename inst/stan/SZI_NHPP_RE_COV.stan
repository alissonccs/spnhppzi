#include /chunks/mylib.stan

data{
  int <lower=0> N;
  int <lower=0>n;
  int <lower=0> p;
  int <lower=0> q;
  // int <lower=0> tp_hf;
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
  matrix [q == 0 ? 0 : N, q] Z;
  matrix [q == 0 ? 0 : n, q] Z1;
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
  real mu_psi;
  real <lower=0> sigma_psi;
  matrix[N,m] g;
  matrix[n,m] G;
  real h1_gamma;
  real h2_gamma;
  real<lower=0> zeta;
    }

parameters{
  // vector  <lower=0> [m] alpha;
  vector <lower=0> [baseline == 2 ? 0 : m]  alpha;
  vector <lower=0> [baseline != 2 ? 0 : m]  gamma;
  vector [p==0 ? 0 :p] beta;
  vector [q==0 ? 0 :q] psi;
  vector[n] omega_1;
  real <lower=0> sigma2_z;
          }

//Estima valor de prob para o modelo que contém somente o intercepto
// transformed parameters{
//   vector [q == 0 ? 0 : q] prob;
//   prob[1]=exp(psi[1])/(1+exp(psi[1]));
//                       }

// transformed parameters {
//  vector [n] log_lik1;
model {
  vector [n] Lambda0 ;
  vector [N] log_lambda0 ;
  vector [N] log_lambda0_event ;
  vector[p == 0 ? 0 : N] eta;
  vector [p == 0 ? 0 : N] eta_event ;
  vector[q == 0 ? 0 : n] eta2;
  //vector [N] eta2_ ;
  vector[n] nu;
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


if(q>0){
     eta2 = Z1*psi;
     }


if(baseline==1){
  Lambda0 = Lambda_plp(max_stop, alpha,n);
  log_lambda0 = log_lambda_plp(time, N, alpha);
  log_lambda0_event = event .*log_lambda0;
}

if(baseline==2){
  Lambda0=Lambda_bp(G, gamma,n);
  log_lambda0=log_lambda_bp(g,gamma, N, zeta);
  log_lambda0_event = event .*log_lambda0;
}


 // CALCULA VEROSSIMILHANÇA ACUMULADA POR INDIVÍDUO
  if(p == 0){
    for (i in 1:n) {

      nu[i] = inv_logit(eta2[i]);

       if(IndRec2[i] == 0)
         target += log_sum_exp(bernoulli_lpmf(1| nu[i]),
                   Lambda0[i] +
                   bernoulli_lpmf(0 | nu[i]));
       else
        target += bernoulli_lpmf(0 | nu[i])+
                   Lambda0[i] +
                   sum(log_lambda0_event[begin_ind[i]:end_ind[i]])
                   // sum_log_lambda0[i]
                   ;
                  }
          }
    else{
      for (i in 1:n) {

        nu[i] = inv_logit(eta2[i]);

       if(IndRec2[i] == 0)
         target += log_sum_exp(bernoulli_lpmf(1| nu[i]),
                   Lambda0[i]*exp_etay[i] +
                   bernoulli_lpmf(0 | nu[i]));
       else
         target +=bernoulli_lpmf(0 | nu[i])+
                  Lambda0[i]*exp_etay[i] +
                  sum(log_lambda0_event[begin_ind[i]:end_ind[i]])+
                  sum(eta_event[begin_ind[i]:end_ind[i]])
                  // sum_log_lambda0[i]+sum_eta[i]
                  ;
                     }
        }

// model{
//      target +=log_lik1;


if (approach==1){
            gamma ~ lognormal(h1_gamma, h2_gamma);
            beta ~ normal(mu_beta,sigma_beta);
            sigma2_z ~ gamma(shp_sigma2_z,scl_sigma2_z);
            omega_1 ~ normal(mu_omega,sigma2_z);
            psi ~ normal(mu_psi,sigma_psi);
  }
}

  generated quantities{
  vector[n] log_lik;
  {

  vector [n] Lambda0 ;
  vector [N] log_lambda0 ;
  vector [N] log_lambda0_event ;
  vector[p == 0 ? 0 : N] eta;
  vector [p == 0 ? 0 : N] eta_event ;
  vector[q == 0 ? 0 : n] eta2;
  //vector [N] eta2_ ;
  vector[n] nu;
  vector[p == 0 ? 0 : n] exp_etay;
  vector [n] sum_log_lambda0 = rep_vector(0, n);
  vector [n] sum_eta = rep_vector(0, n);
  int a = 0;
  int c = 0;

 if(p>0 ){
  for (i in 1:N){
     eta[i] = X[i,]*beta+omega_1[id[i]];
     eta_event[i] = event[i]*eta[i];
  }
  for (j in 1:n){
     exp_etay[j] = exp(Xy[j,]*beta+omega_1[j]);
        }
}

if(q>0){
     eta2 = Z1*psi;
     }


if(baseline==1){
  Lambda0 = Lambda_plp(max_stop, alpha,n);
  log_lambda0 = log_lambda_plp(time, N, alpha);
  log_lambda0_event = event .*log_lambda0;
}


if(baseline==2){
  Lambda0=Lambda_bp(G, gamma,n);
  log_lambda0=log_lambda_bp(g,gamma, N, zeta);
  log_lambda0_event = event .*log_lambda0;
}


 // CALCULA VEROSSIMILHANÇA ACUMULADA POR INDIVÍDUO
 // for (b in 1:n) {
 //    sum_log_lambda0[b] = 0;
 //    sum_eta[b] = 0;
 //  for (i in 1:n_ind[b]) {
 //    if (a == N) {
 //        break;
 //                }
 //    a = a + 1;
 //    sum_log_lambda0[b] += log_lambda0_event[a];
 //    if(p>0){
 //    sum_eta[b] += eta_event[a];
 //    }
 //                         }
 //    if (a == N) {
 //        break;
 //    }
 //  }

//
// for ( b in 1:n) {
//         sum_log_lambda0[b]=sum(log_lambda0_event[begin_ind[b]:end_ind[b]]);
//         if(p>0){
//        sum_eta[b]=sum(eta_event[begin_ind[b]:end_ind[b]]);
//                }
//  }


 if(p == 0){
    for (i in 1:n) {

      nu[i] = inv_logit(eta2[i]);

       if(IndRec2[i] == 0)
        log_lik[i]= log_sum_exp(bernoulli_lpmf(1| nu[i]),
                   Lambda0[i] +
                   bernoulli_lpmf(0 | nu[i]));
       else
        log_lik[i]= bernoulli_lpmf(0 | nu[i])+
                   Lambda0[i] +
                   sum(log_lambda0_event[begin_ind[i]:end_ind[i]])
                   // sum_log_lambda0[i]
                   ;
                  }
          }
    else{
      for (i in 1:n) {

        nu[i] = inv_logit(eta2[i]);

       if(IndRec2[i] == 0)
       log_lik[i]= log_sum_exp(bernoulli_lpmf(1| nu[i]),
                   Lambda0[i]*exp_etay[i] +
                   bernoulli_lpmf(0 | nu[i]));
       else
       log_lik[i]=bernoulli_lpmf(0 | nu[i])+
                  Lambda0[i]*exp_etay[i] +
                  sum(log_lambda0_event[begin_ind[i]:end_ind[i]])+
                  sum(eta_event[begin_ind[i]:end_ind[i]])
                  // sum_log_lambda0[i]+sum_eta[i]
                  ;
                     }
        }
  }
  }
