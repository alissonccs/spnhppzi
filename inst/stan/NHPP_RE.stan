#include /chunks/mylib.stan
data{
  int <lower=0> N;
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
    }

parameters{
  vector  <lower=0> [m] alpha;
  vector [p] beta;
  real <lower=0> sigma2_z;
  vector [n] omega;
            }

transformed parameters {
  vector[n] log_lik;
  {

  vector [n] Lambda0 ;
  vector [N] log_lambda0 ;
  vector [N] log_lambda0_event;
  vector[p == 0 ? 0 : N] eta;
  vector [N] eta_event ;
  vector [p == 0 ? 0 : n] exp_etay;
  vector [n] sum_log_lambda0 = rep_vector(0, n);
  vector [n] sum_eta = rep_vector(0, n);
  int a = 0;
  int c = 0;

if(p>0){
  for (i in 1:N){
     eta[i] = X[i,]*beta+omega[id[i]];
     eta_event[i] = event[i]*eta[i];
  }
     exp_etay = exp(Xy*beta+omega);
        }

  Lambda0 = Lambda_plp(max_stop, alpha,n);
  log_lambda0 = log_lambda_plp(time, N, alpha);
  log_lambda0_event= event .*log_lambda0;

 for ( b in 1:n) {
        sum_log_lambda0[b]=sum(log_lambda0_event[begin_ind[b]:end_ind[b]]);
        if(p>0){
       sum_eta[b]=sum(eta_event[begin_ind[b]:end_ind[b]]);
               }
 }


if(p == 0){
    for (i in 1:n) {
      log_lik[i]=   Lambda0[i] +
                   sum_log_lambda0[i];
                  }
          }
    else{
      for (i in 1:n) {
      log_lik[i]= Lambda0[i]*exp_etay[i]+
                 sum_log_lambda0[i]+sum_eta[i];

                     }
        }


  //   if(p==0){
  //    log_lik=sum(Lambda0) +
  //    sum(event .*log_lambda0);
  //         }
  // else{
  //    log_lik= sum(Lambda0 .* exp_etay) +
  //    sum(event .*(log_lambda0 + eta));
  //      }
   }
}
model{
     target +=log_lik;

if(approach==1){
            alpha[1] ~ gamma(shp_alpha1,scl_alpha1);
            alpha[2] ~ gamma(shp_alpha2,scl_alpha2);
            beta ~ normal(mu_beta,sigma_beta);
            sigma2_z ~ gamma(shp_sigma2_z,scl_sigma2_z);
            omega~ normal(log(1 / sqrt(sigma2_z + 1)),sqrt(log(sigma2_z + 1)));
                                         }
  }
