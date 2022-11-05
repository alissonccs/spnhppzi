#include /chunks/mylib.stan

data{
  int <lower=0> N;
  int <lower=0>n;
  int <lower=0> p;
  int <lower=0> q;
  int <lower=0> tp_hf;
  int m;
  int n_ind [n];
  int begin_int [n];
  int end_int[n];
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
 vector [n] log_lik1;
 {
 vector [n] Lambda0 ;
 vector [N] log_lambda0 ;
 vector [N] log_lambda0_ ;
 vector[p == 0 ? 0 : N] eta1;
 vector [p == 0 ? 0 : N] eta1_ ;
 vector[q == 0 ? 0 : n] eta2;
 //vector [N] eta2_ ;
 vector[n] omega;
 vector[p == 0 ? 0 : n] exp_etay;
 vector [n] sum_log_lambda0 = rep_vector(0, n);
 vector [n] sum_eta1 = rep_vector(0, n);
 int a = 0;

if(p>0){
    eta1 = X*beta;
    eta1_ = event .*eta1;
    exp_etay = exp(Xy*beta);
       }
if(q>0){
    eta2 = Z1*psi;
    }
 Lambda0 = Lambda_plp2(max_stop, alpha,n);
 log_lambda0 = log_lambda_plp2(time, N, alpha);
 log_lambda0_ = event .*log_lambda0;
// CALCULA VEROSSIMILHANÇA ACUMULADA POR INDIVÍDUO
// for (b in 1:n) {
//    sum_log_lambda0[b] = 0;
//    sum_eta[b] = 0;
//  for (i in 1:n_ind[b]) {
//    if (a == N) {
//        break;
//                }
//    a = a + 1;
//    sum_log_lambda0[b] += log_lambda0_[a];
//    if(p>0){
//    sum_eta[b] += eta1_[a];
//    }
//                         }
//    if (a == N) {
//        break;
//    }
//  }
for ( b in 1:n) {
       sum_log_lambda0[b]=sum(log_lambda0_[begin_int[b]:end_int[b]]);
       if(p>0){
      sum_eta1[b]=sum(eta1_[begin_int[b]:end_int[b]]);
              }
}
if(p == 0){
   for (i in 1:n) {
     omega[i] = inv_logit(eta2[i]);
      if(IndRec2[i] == 0)
        log_lik1[i]= log_sum_exp(bernoulli_lpmf(1| omega[i]),
                  Lambda0[i] +
                  bernoulli_lpmf(0 | omega[i]));
      else
        log_lik1[i]= bernoulli_lpmf(0 | omega[i])+
                  Lambda0[i] +
                  sum_log_lambda0[i];
                 }
         }
   else{
     for (i in 1:n) {
       omega[i] = inv_logit(eta2[i]);
      if(IndRec2[i] == 0)
        log_lik1[i]= log_sum_exp(bernoulli_lpmf(1| omega[i]),
                  Lambda0[i]*exp_etay[i] +
                  bernoulli_lpmf(0 | omega[i]));
      else
        log_lik1[i]+=bernoulli_lpmf(0 | omega[i])+
                 Lambda0[i]*exp_etay[i] +
                 sum_log_lambda0[i]+sum_eta1[i];
                    }
       }
 }
}


// model{
//   vector [n] Lambda0 ;
//   vector [N] log_lambda0 ;
//   vector [N] log_lambda0_ ;
//   vector[p == 0 ? 0 : N] eta1;
//   vector [p == 0 ? 0 : N] eta1_ ;
//   vector[q == 0 ? 0 : n] eta2;
//   //vector [N] eta2_ ;
//   vector[n] omega;
//   vector[p == 0 ? 0 : n] exp_etay;
//   vector [n] sum_log_lambda0 = rep_vector(0, n);
//   vector [n] sum_eta1 = rep_vector(0, n);
//   int a = 0;
//
//  if(p>0){
//      eta1 = X*beta;
//      eta1_ = event .*eta1;
//      exp_etay = exp(Xy*beta);
//         }
//
//  if(q>0){
//      eta2 = Z1*psi;
//      }
//
//
//   Lambda0 = Lambda_plp2(max_stop, alpha,n);
//   log_lambda0 = lambda_plp2(time, N, alpha);
//   log_lambda0_ = event .*log_lambda0;


 // CALCULA VEROSSIMILHANÇA ACUMULADA POR INDIVÍDUO
 // for (b in 1:n) {
 //    sum_log_lambda0[b] = 0;
 //    sum_eta[b] = 0;
 //  for (i in 1:n_ind[b]) {
 //    if (a == N) {
 //        break;
 //                }
 //    a = a + 1;
 //    sum_log_lambda0[b] += log_lambda0_[a];
 //    if(p>0){
 //    sum_eta[b] += eta1_[a];
 //    }
 //                         }
 //    if (a == N) {
 //        break;
 //    }
 //  }

 // for ( b in 1:n) {
 //        sum_log_lambda0[b]=sum(log_lambda0_[begin_int[b]:end_int[b]]);
 //        if(p>0){
 //       sum_eta1[b]=sum(eta1_[begin_int[b]:end_int[b]]);
 //               }
 // }
 //
 // if(p == 0){
 //    for (i in 1:n) {
 //
 //      omega[i] = inv_logit(eta2[i]);
 //
 //       if(IndRec2[i] == 0)
 //         target += log_sum_exp(bernoulli_lpmf(1| omega[i]),
 //                   Lambda0[i] +
 //                   bernoulli_lpmf(0 | omega[i]));
 //       else
 //         target += bernoulli_lpmf(0 | omega[i])+
 //                   Lambda0[i] +
 //                   sum_log_lambda0[i];
 //                  }
 //          }
 //    else{
 //      for (i in 1:n) {
 //
 //        omega[i] = inv_logit(eta2[i]);
 //
 //       if(IndRec2[i] == 0)
 //         target += log_sum_exp(bernoulli_lpmf(1| omega[i]),
 //                   Lambda0[i]*exp_etay[i] +
 //                   bernoulli_lpmf(0 | omega[i]));
 //       else
 //         target +=bernoulli_lpmf(0 | omega[i])+
 //                  Lambda0[i]*exp_etay[i] +
 //                  sum_log_lambda0[i]+sum_eta1[i];
 //                     }
 //        }
 //
 // if(approach==1 && tp_prior==1){
 //            alpha[1] ~ gamma(shp_alpha1,scl_alpha1);
 //            alpha[2] ~ gamma(shp_alpha2,scl_alpha2);
 //            beta ~ normal(mu_beta,sigma_beta);
 //            psi ~ normal(mu_psi,sigma_psi);
 //                                          }
 //  }

model{
     target +=log_lik1;

if(approach==1 && tp_prior==1){
            alpha[1] ~ gamma(shp_alpha1,scl_alpha1);
            alpha[2] ~ gamma(shp_alpha2,scl_alpha2);
            beta ~ normal(mu_beta,sigma_beta);
            psi ~ normal(mu_psi,sigma_psi);
                              }
  }
