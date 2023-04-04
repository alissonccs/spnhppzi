#include /chunks/mylib.stan

data{
  int <lower=0> N;
  int <lower=0> p;
  int <lower=0> n;
  int <lower=0> baseline;
  int m;
  int n_ind [n];
  int id[N];
  // int begin_int [n];
  // int end_int[n];
  vector [N] event;
  vector [n] max_stop;
  vector [N] time;
  matrix [p == 0 ? 0 : N, p] X;
  matrix [p == 0 ? 0 : n, p] Xy;
  int<lower=0> approach;
  real shp_alpha1;
  real scl_alpha1;
  real shp_alpha2;
  real scl_alpha2;
  real mu_beta;
  real <lower=0> sigma_beta;
  int <lower=0> tp_prior;
 }

parameters{
  vector  <lower=0> [m] alpha;
  vector [p] beta;
  }

model{vector [n] Lambda0 ;
  vector [N] log_lambda0 ;
  vector[p == 0 ? 0 : N] eta;
  vector[p == 0 ? 0 : n] exp_etay;

 if(p>0){
     eta = X*beta;
     exp_etay = exp(Xy*beta);
         }

  Lambda0=Lambda_plp2(max_stop, alpha,n);
  log_lambda0=log_lambda_plp2(time, N, alpha);



  if(p==0){
     target +=sum(Lambda0) +
     sum(event .*log_lambda0);
          }
  else{
     target += sum(Lambda0 .* exp_etay) +
     sum(event .*(log_lambda0 + eta));
       }

 if(approach==1 && tp_prior==1){
     alpha[1] ~ gamma(shp_alpha1,scl_alpha1);
     alpha[2] ~ gamma(shp_alpha2,scl_alpha2);
     beta ~ normal(mu_beta,sigma_beta);
                               }
  }

  generated quantities {
  row_vector [n] log_lik;
  {
  vector [n] Lambda0 ;
  vector [N] log_lambda0;
  vector [N] log_lambda0_event;
  vector[p == 0 ? 0 : N] eta;
  vector [p == 0 ? 0 :N] eta_event;
  vector[p == 0 ? 0 : n] exp_etay;
  vector [n] sum_log_lambda0 = rep_vector(0, n);
  vector [n] sum_eta = rep_vector(0, n);
  int a = 0;
  // int c = 0;

  if(p>0){
  for (i in 1:N){
     eta[i] = X[i,]*beta;
     eta_event[i] = event[i]*eta[i];
  }
     exp_etay = exp(Xy*beta);
        }

  Lambda0 = Lambda_plp2(max_stop, alpha,n);
  log_lambda0 = log_lambda_plp2(time, N, alpha);
  log_lambda0_event = event .*log_lambda0;

   // CALCULA VEROSSIMILHANÇA ACUMULADA POR INDIVÍDUO
    for (b in 1:n) {
       sum_log_lambda0[b] = 0;
       sum_eta[b] = 0;
     for (i in 1:n_ind[b]) {
       if (a <= N) {
       a = a + 1;
       sum_log_lambda0[b] += log_lambda0_event[a];
       if(p>0){
       sum_eta[b] += eta_event[a];
       }
                    }
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
  }
}
