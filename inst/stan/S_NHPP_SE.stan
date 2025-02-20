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
  real shp_sigma2_z;
  real scl_sigma2_z;
  real mu_beta;
  real <lower=0> sigma_beta;
  int SP_ID[N];
  int gr_SP_ID[n];
  matrix<lower = 0, upper = 1>[SP_N,SP_N] nb_mat; // adjacency matrix
  int W_n;                // number of adjacent region pairs
  real shp_tau;
  real scl_tau;
  matrix[N,m] g;
  matrix[n,m] G;
  real h1_gamma;
  real h2_gamma;
  real<lower=0> zeta;
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
  vector [p] beta;
  vector <lower=0> [baseline == 2 ? 0 : m]  alpha;
  vector <lower=0> [baseline != 2 ? 0 : m]  gamma;
    vector [SP_N] omega;
  real<lower = 0> tau;

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
  Lambda0 = Lambda_plp(max_stop, alpha,n);
  log_lambda0 = log_lambda_plp(time, N, alpha);
  log_lambda0_event = event .*log_lambda0;
}

if(baseline==2){
  Lambda0=Lambda_bp(G, gamma,n);
  log_lambda0=log_lambda_bp(g,gamma, N, zeta);
  log_lambda0_event = event .*log_lambda0;
}

 //  for ( b in 1:n) {
 //        sum_log_lambda0[b]=sum(log_lambda0_event[begin_ind[b]:end_ind[b]]);
 //        if(p>0){
 //       sum_eta[b]=sum(eta_event[begin_ind[b]:end_ind[b]]);
 //               }
 // }


 if(p == 0){ for (i in 1:n) {
      target +=   Lambda0[i] +
                   sum(log_lambda0_event[begin_ind[i]:end_ind[i]]);
                  }
          }
    else{
      for (i in 1:n) {
      target += Lambda0[i]*exp_etay[i]+
                 sum(log_lambda0_event[begin_ind[i]:end_ind[i]])+
                 sum(eta_event[begin_ind[i]:end_ind[i]]);

                     }
        }

if(approach==1 && baseline==2){
            gamma ~ lognormal(h1_gamma, h2_gamma);
            beta ~ normal(mu_beta,sigma_beta);
            omega ~ sparse_iar(tau, W_sparse, D_sparse, lambda, SP_N, W_n);
            sum(omega) ~ normal(0, 0.001 * SP_N);
            tau ~ gamma(shp_tau, scl_tau);
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
     eta[i] = X[i,]*beta+omega[SP_ID[i]];
     eta_event[i] = event[i]*eta[i];
  }
     for (j in 1:n){
      exp_etay[j] = exp(Xy[j,]*beta+omega[gr_SP_ID[j]]);
    }
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

 //  for ( b in 1:n) {
 //        sum_log_lambda0[b]=sum(log_lambda0_event[begin_ind[b]:end_ind[b]]);
 //        if(p>0){
 //       sum_eta[b]=sum(eta_event[begin_ind[b]:end_ind[b]]);
 //               }
 // }


 if(p == 0){
   for (i in 1:n) {
      log_lik[i]=   Lambda0[i] +
                   sum(log_lambda0_event[begin_ind[i]:end_ind[i]]);
                  }
          }
    else{
      for (i in 1:n) {
      log_lik[i]= Lambda0[i]*exp_etay[i]+
                 sum(log_lambda0_event[begin_ind[i]:end_ind[i]])+
                 sum(eta_event[begin_ind[i]:end_ind[i]]);

                     }
        }
  }
}
