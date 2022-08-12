#include /chunks/mylib.stan
data{
  int <lower=0> N;
  int <lower=0> p;
  int <lower=0>n;
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
  // matrix[n,3] position_ind;
  vector [n] IndRec2;
  int<lower=0> approach;
  int <lower=0> ZI;
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
  real <lower=0,upper=1> pii [ZI == 0 ? 0 : 1];
          }

model{
  vector [n] Lambda0 ;
  vector [N] log_lambda0 ;
  vector [N] log_lambda0_event ;
  vector[p == 0 ? 0 : N] eta;
  vector [N] eta_event ;
  vector[p == 0 ? 0 : n] exp_etay;
  vector [n] sum_log_lambda0 = rep_vector(0, n);
  vector [n] sum_eta = rep_vector(0, n);
  int a = 0;

 if(p>0){
     eta = X*beta;
     eta_event = event .*eta;
     exp_etay = exp(Xy*beta);
        }

  Lambda0 = Lambda_plp2(max_stop, alpha,n);
  log_lambda0 = lambda_plp2(time, N, alpha);
  log_lambda0_event = event .*log_lambda0;


 for ( b in 1:n) {
        sum_log_lambda0[b]=sum(log_lambda0_event[begin_ind[b]:end_ind[b]]);
        if(p>0){
       sum_eta[b]=sum(eta_event[begin_ind[b]:end_ind[b]]);
               }
 }

 if(p == 0){
    // for (i in 1:n) {

       // if(IndRec2[i] == 0)
         target +=  (1-IndRec2).* (log(exp(bernoulli_lpmf(1| pii))+
                   exp(Lambda0 +
                   bernoulli_lpmf(0 | pii))));
       // else
         target += IndRec2 .*(bernoulli_lpmf(0 | pii)+
                   Lambda0 +
                   sum_log_lambda0);
                  // }
          }
    else{
      // for (i in 1:n) {
      //  if(IndRec2[i] == 0)
         target +=sum((1-IndRec2).* (log(exp(bernoulli_lpmf(1| pii))+
                   exp(Lambda0 .* exp_etay +
                   bernoulli_lpmf(0 | pii)))));
       // else
       // for (j in 1:N[id[i]]){
         target +=  sum(IndRec2 .*(bernoulli_lpmf(0 | pii)+
                    Lambda0 .* exp_etay +
                    sum_log_lambda0+
                    sum_eta));
                     }
       // }

 if(approach==1 && tp_prior==1){
            alpha[1] ~ gamma(shp_alpha1,scl_alpha1);
            alpha[2] ~ gamma(shp_alpha2,scl_alpha2);
            beta ~ normal(mu_beta,sigma_beta);
                               }
  }
