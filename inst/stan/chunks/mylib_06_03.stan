functions{

vector Lambda_plp1(vector max_stop, vector alpha, int n){
      vector [n] loglik;
      for(i in 1:n){
     loglik[i] =  -exp(lmultiply(alpha[2],max_stop[i]) - lmultiply(alpha[2],alpha[1]));
           }
      return loglik;
    }
vector log_lambda_plp1(vector time, int N, vector alpha){
      vector[N] loglik1;
      for(i in 1:N){
      loglik1[i] = log(alpha[2])-
      lmultiply(alpha[2], alpha[1])+
      lmultiply(alpha[2]-1,time[i]);
      // lprob1[i] = log(alpha[1])+log(alpha[2])+
      // lmultiply(alpha[2]-1,time[i]);
      }
      return loglik1;
    }

vector Lambda_plp2(vector max_stop, vector alpha, int n){
      vector [n] loglik;
      for(i in 1:n){
     // lprob[i] =  -exp(lmultiply(alpha[2],max_stop[i]) - lmultiply(alpha[2],alpha[1]));
      loglik[i]= -exp(log(alpha[1])+lmultiply(alpha[2],max_stop[i]));
      }
      return loglik;
    }
vector log_lambda_plp2(vector time, int N, vector alpha){
      vector[N] loglik1;
      for(i in 1:N){
      // lprob1[i] = log(alpha[2])-
      // lmultiply(alpha[2], alpha[1])+
      // lmultiply(alpha[2]-1,time[i]);
      loglik1[i] = log(alpha[1])+log(alpha[2])+
      lmultiply(alpha[2]-1,time[i]);
      }
      return loglik1;
    }

vector Lambda_plp3(vector max_stop, vector alpha, int n, real zeta){
      vector [n] loglik;
      vector [n] y1=max_stop/zeta;
    for(i in 1:n){
     loglik[i]= -exp(log(alpha[1])+lmultiply(alpha[2],y1[i]));
     // lprob[i] =  -exp(lmultiply(alpha[2],y1[i]) - lmultiply(alpha[2],alpha[1]));
           }
      return loglik;
    }
vector log_lambda_plp3(vector time, int N, vector alpha, real zeta){
      vector[N] loglik1;
      vector [N] y=time/zeta;
      for(i in 1:N){
      // lprob1[i] = log(alpha[2])-
      // lmultiply(alpha[2], alpha[1])+
      // lmultiply(alpha[2]-1,y[i]);
      loglik1[i] = log(alpha[1])+log(alpha[2])+
      lmultiply(alpha[2]-1,y[i]);
      }
      return (loglik1-log(zeta));
    }

real sparse_car_lpdf(vector omega, real sp_tau, real sp_alpha,
    int[,] W_sparse, vector D_sparse, vector sp_lambda, int SP_N, int W_n) {
      row_vector[SP_N] omegat_D; // omega' * D
      row_vector[SP_N] omegat_W; // omega' * W
      vector[SP_N] ldet_terms;

      omegat_D = (omega .* D_sparse)';
      omegat_W = rep_row_vector(0, SP_N);
      for (i in 1:W_n) {
        omegat_W[W_sparse[i, 1]] = omegat_W[W_sparse[i, 1]] + omega[W_sparse[i, 2]];
        omegat_W[W_sparse[i, 2]] = omegat_W[W_sparse[i, 2]] + omega[W_sparse[i, 1]];
      }

      for (i in 1:SP_N) ldet_terms[i] = log1m(sp_alpha * sp_lambda[i]);
      return 0.5 * (SP_N * log(sp_tau)
                    + sum(ldet_terms)
                    - sp_tau * (omegat_D * omega - sp_alpha * (omegat_W * omega)));
  }

real sparse_iar_lpdf(vector omega, real sp_tau,
    int[,] W_sparse, vector D_sparse, vector sp_lambda, int SP_N, int W_n) {
      row_vector[SP_N] omegat_D; // omega' * D
      row_vector[SP_N] omegat_W; // omega' * W
      vector[SP_N] ldet_terms;

      omegat_D = (omega .* D_sparse)';
      omegat_W = rep_row_vector(0, SP_N);
      for (i in 1:W_n) {
        omegat_W[W_sparse[i, 1]] = omegat_W[W_sparse[i, 1]] + omega[W_sparse[i, 2]];
        omegat_W[W_sparse[i, 2]] = omegat_W[W_sparse[i, 2]] + omega[W_sparse[i, 1]];
      }

      return 0.5 * ((SP_N-1) * log(sp_tau)
                    - sp_tau * (omegat_D * omega - (omegat_W * omega)));
  }

real sparse_iar1_lpdf(vector omega, real sp_std,
    int[,] W_sparse, vector D_sparse, vector sp_lambda, int SP_N, int W_n) {
      row_vector[SP_N] omegat_D; // omega' * D
      row_vector[SP_N] omegat_W; // omega' * W
      vector[SP_N] ldet_terms;

      omegat_D = (omega .* D_sparse)';
      omegat_W = rep_row_vector(0, SP_N);
      for (i in 1:W_n) {
        omegat_W[W_sparse[i, 1]] = omegat_W[W_sparse[i, 1]] + omega[W_sparse[i, 2]];
        omegat_W[W_sparse[i, 2]] = omegat_W[W_sparse[i, 2]] + omega[W_sparse[i, 1]];
      }

      return 0.5 * (-(SP_N-1) * 2*log(sp_std)
                    - (1/sp_std^2) * (omegat_D * omega - (omegat_W * omega)));
  }

vector Lambda_bp(matrix G, vector gamma,int n){
      vector [n] loglik;
      loglik= -G*gamma;
      return loglik;
    }
vector log_lambda_bp(matrix g, vector gamma,int N,real zeta){
      vector[N] loglik1;
      loglik1 = log(g*gamma) - log(zeta);
            return loglik1;
    }
}

