#include /chunks/mylib.stan
data{
  int <lower=0> N;
  int <lower=0>n;
  int <lower=0> p;
  int <lower=0> q;
  int <lower=0> tp_hf;
  int m;
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
  real mu_beta;
  real <lower=0> sigma_beta;
  real mu_psi;
  real <lower=0> sigma_psi;
  int <lower=0> tp_prior;
    }

parameters{
  vector  <lower=0> [m] alpha;
  vector [p==0 ? 0 :p] beta;
  vector [q==0 ? 0 :q] psi;
  //real <lower=0,upper=1> prob [q == 0 ? 1 : 0];
          }

//Estima valor de prob para o modelo que contém somente o intercepto
// transformed parameters{
//   vector [q == 0 ? 0 : q] prob;
//   prob[1]=exp(psi[1])/(1+exp(psi[1]));
//                       }

transformed parameters {
 vector [n] log_lik;
 {
 // model{
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

if(p>0){
    eta = X*beta;
    eta_event = event .*eta;
    exp_etay = exp(Xy*beta);
       }
if(q>0){
    eta2 = Z1*psi;
    }
 Lambda0 = Lambda_plp2(max_stop, alpha,n);
 log_lambda0 = log_lambda_plp2(time, N, alpha);
 log_lambda0_event = event .*log_lambda0;
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
// for ( b in 1:n) {
//        sum_log_lambda0[b]=sum(log_lambda0_event[begin_ind[b]:end_ind[b]]);
//        if(p>0){
//       sum_eta[b]=sum(eta_event[begin_ind[b]:end_ind[b]]);
//               }
// }
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
      log_lik[i]= bernoulli_lpmf(0 | nu[i])+
                 Lambda0[i]*exp_etay[i] +
                 sum(log_lambda0_event[begin_ind[i]:end_ind[i]])+
                 sum(eta_event[begin_ind[i]:end_ind[i]])
                 // sum_log_lambda0[i]+sum_eta[i]
                 ;
                    }
       }
 }
}
//
model{
     target +=log_lik;

if(approach==1 && tp_prior==1){
            alpha[1] ~ gamma(shp_alpha1,scl_alpha1);
            alpha[2] ~ gamma(shp_alpha2,scl_alpha2);
            beta ~ normal(mu_beta,sigma_beta);
            psi ~ normal(mu_psi,sigma_psi);
                              }
  }
