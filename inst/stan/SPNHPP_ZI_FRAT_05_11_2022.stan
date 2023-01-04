#include /chunks/spmylib.stan
data{
  int <lower=0> N;
  int <lower=0> SP_N;
  int <lower=0> p;
  int <lower=0>n;
  // int <lower=0> tp_hf;
  int <lower=0> baseline;
  int m;
  int id[N];
  int n_ind [n];
  int begin_ind [n];
  int end_ind[n];
  vector [n] n_ind1;
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
  // real shp_sigma_omega;
  // real scl_sigma_omega;
  real shp_sigma2_z;
  real scl_sigma2_z;
  real mu_beta;
  real <lower=0> sigma_beta;
  int <lower=0> tp_prior;
  int SP_ID[N];
  int gr_SP_ID[n];
  matrix<lower = 0, upper = 1>[SP_N,SP_N] nb_mat; // adjacency matrix
  int W_n;                // number of adjacent region pairs
  real shp_tau;
  real scl_tau;
    }

transformed data{

  int W_sparse[W_n, 2];   // adjacency pairs
  vector[SP_N] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[SP_N] lambda;       // eigenvalues of invsqrtD * W * invsqrtD

  { // generate sparse representation for W
  int counter;
  counter = 1;
  // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1:(SP_N - 1)) {
      for (j in (i + 1):SP_N) {
        if (nb_mat[i, j] == 1) {
          W_sparse[counter, 1] = i;
          W_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  for (i in 1:SP_N) D_sparse[i] = sum(nb_mat[i]);
  {
    vector[SP_N] invsqrtD;
    for (i in 1:SP_N) {
      invsqrtD[i] = 1 / sqrt(D_sparse[i]);
    }
    lambda = eigenvalues_sym(quad_form(nb_mat, diag_matrix(invsqrtD)));
  }
}

parameters{
  vector  <lower=0> [m] alpha;
  vector [p] beta;
  real <lower=0,upper=1> pii [ZI == 0 ? 0 : 1];
  vector [SP_N] omega;
  real<lower = 0> tau;
  // vector [n] omega;
  // real <lower=0> sigma_omega;
  // real <lower=0> sigma2_z;
          }


// transformed parameters {
  // vector[n] log_lik1;
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
     eta[i] = X[i,]*beta+omega[SP_ID[i]];
     eta_event[i] = event[i]*eta[i];
  }
     for (j in 1:n){
      exp_etay[j] = exp(Xy[j,]*beta+omega[gr_SP_ID[j]]);
    }
        }


if(baseline==1){
  Lambda0 = Lambda_plp1(max_stop, alpha,n);
  log_lambda0 = log_lambda_plp1(time, N, alpha);
  log_lambda0_event = event .*log_lambda0;
}

if(baseline==2){
  Lambda0 = Lambda_plp2(max_stop, alpha,n);
  log_lambda0 = log_lambda_plp2(time, N, alpha);
  log_lambda0_event = event .*log_lambda0;
}

if(baseline==3){
  Lambda0 = Lambda_plp3(max_stop, alpha,n,tau);
  log_lambda0 = log_lambda_plp3(time, N, alpha,tau);
  log_lambda0_event = event .*log_lambda0;
}

// if(baseline==4){
//   Lambda0=Lambda_bp(G, gamma,n);
//   log_lambda0=log_lambda_bp(g,gamma, N, tau);
//   log_lambda0_event = event .*log_lambda0;
// }

 //  for ( b in 1:n) {
 //        sum_log_lambda0[b]=sum(log_lambda0_event[begin_ind[b]:end_ind[b]]);
 //        if(p>0){
 //       sum_eta[b]=sum(eta_event[begin_ind[b]:end_ind[b]]);
 //               }
 // }


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

 // if(approach==1 && tp_prior==1){
 //            alpha[1] ~ gamma(shp_alpha1,scl_alpha1);
 //            alpha[2] ~ gamma(shp_alpha2,scl_alpha2);
 //            beta ~ normal(mu_beta,sigma_beta);
 //            // sigma_omega ~ gamma(shp_sigma_omega,scl_sigma_omega);
 //            omega[id] ~ normal(mu_omega,sigma_omega);
 //                               }
//   }
// }
// model{
//      target +=log_lik1;

if(approach==1 && tp_prior==1 && (baseline==1 || baseline==2)){
            alpha[1] ~ gamma(shp_alpha1,scl_alpha1);
            alpha[2] ~ gamma(shp_alpha2,scl_alpha2);
            beta ~ normal(mu_beta,sigma_beta);
            // sigma2_z ~ gamma(shp_sigma2_z,scl_sigma2_z);
            // sigma_omega ~ gamma(shp_sigma2_z,scl_sigma2_z);
            // omega~ normal(-(sigma_omega)^2/2,sigma_omega);
            // sigma_omega ~ gamma(shp_sigma_omega,scl_sigma_omega);
            // omega~ normal(log(1 / sqrt(sigma2_z + 1)),sqrt(log(sigma2_z + 1)));
            // omega ~ normal(mu_omega,sigma_omega);
            omega ~sparse_iar_lpdf(tau, W_sparse, D_sparse, lambda, SP_N, W_n);
            tau ~ gamma(shp_tau, scl_tau);
            sum(omega) ~ normal(0, 0.001 * SP_N);
                               }

  if(approach==1 && tp_prior==1 && baseline==3){
            alpha[1] ~ lognormal(shp_alpha1,scl_alpha1);
            alpha[2] ~ gamma(shp_alpha2,scl_alpha2);
            beta ~ normal(mu_beta,sigma_beta);
            // sigma2_z ~ gamma(shp_sigma2_z,scl_sigma2_z);
            // sigma_omega ~ gamma(shp_sigma2_z,scl_sigma2_z);
            // omega~ normal(-(sigma_omega)^2/2,sigma_omega);
            // sigma_omega ~ gamma(shp_sigma_omega,scl_sigma_omega);
            // omega~ normal(log(1 / sqrt(sigma2_z + 1)),sqrt(log(sigma2_z + 1)));
            // omega ~ normal(mu_omega,sigma_omega);
            omega ~sparse_iar_lpdf(tau, W_sparse, D_sparse, lambda, SP_N, W_n);
            tau ~ gamma(shp_tau, scl_tau);
            sum(omega) ~ normal(0, 0.001 * SP_N);
                               }
  }
