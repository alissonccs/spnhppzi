functions{

vector Lambda_plp(vector max_stop, vector alpha, int n){
      vector [n] lprob;
      for(i in 1:n){
     // lprob[i] =  -exp(lmultiply(alpha[2],max_stop[i]) - lmultiply(alpha[2],alpha[1]));
      lprob[i]= -exp(log(alpha[1])+lmultiply(alpha[2],max_stop[i]));
      }
      return lprob;
    }
vector log_lambda_plp(vector time, int N, vector alpha){
      vector[N] lprob1;
      for(i in 1:N){
      // lprob1[i] = log(alpha[2])-
      // lmultiply(alpha[2], alpha[1])+
      // lmultiply(alpha[2]-1,time[i]);
      lprob1[i] = log(alpha[1])+log(alpha[2])+
      lmultiply(alpha[2]-1,time[i]);
      }
      return lprob1;
    }


real sparse_car_lpdf(vector omega, real tau, real alpha,
    int[,] W_sparse, vector D_sparse, vector lambda, int SP_N, int W_n) {
      row_vector[SP_N] omegat_D; // omega' * D
      row_vector[SP_N] omegat_W; // omega' * W
      vector[SP_N] ldet_terms;

      omegat_D = (omega .* D_sparse)';
      omegat_W = rep_row_vector(0, SP_N);
      for (i in 1:W_n) {
        omegat_W[W_sparse[i, 1]] = omegat_W[W_sparse[i, 1]] + omega[W_sparse[i, 2]];
        omegat_W[W_sparse[i, 2]] = omegat_W[W_sparse[i, 2]] + omega[W_sparse[i, 1]];
      }

      for (i in 1:SP_N) ldet_terms[i] = log1m(alpha * lambda[i]);
      return 0.5 * (SP_N * log(tau)
                    + sum(ldet_terms)
                    - tau * (omegat_D * omega - alpha * (omegat_W * omega)));
  }

real sparse_iar_lpdf(vector omega, real tau,
    int[,] W_sparse, vector D_sparse, vector lambda, int SP_N, int W_n) {
      row_vector[SP_N] omegat_D; // omega' * D
      row_vector[SP_N] omegat_W; // omega' * W
      vector[SP_N] ldet_terms;

      omegat_D = (omega .* D_sparse)';
      omegat_W = rep_row_vector(0, SP_N);
      for (i in 1:W_n) {
        omegat_W[W_sparse[i, 1]] = omegat_W[W_sparse[i, 1]] + omega[W_sparse[i, 2]];
        omegat_W[W_sparse[i, 2]] = omegat_W[W_sparse[i, 2]] + omega[W_sparse[i, 1]];
      }

      return 0.5 * ((SP_N-1) * log(tau)
                    - tau * (omegat_D * omega - (omegat_W * omega)));
  }

real sparse_iar1_lpdf(vector omega, real tau,
    int[,] W_sparse, vector D_sparse, vector lambda, int SP_N, int W_n) {
      row_vector[SP_N] omegat_D; // omega' * D
      row_vector[SP_N] omegat_W; // omega' * W
      vector[SP_N] ldet_terms;

      omegat_D = (omega .* D_sparse)';
      omegat_W = rep_row_vector(0, SP_N);
      for (i in 1:W_n) {
        omegat_W[W_sparse[i, 1]] = omegat_W[W_sparse[i, 1]] + omega[W_sparse[i, 2]];
        omegat_W[W_sparse[i, 2]] = omegat_W[W_sparse[i, 2]] + omega[W_sparse[i, 1]];
      }

      return 0.5 * (-(SP_N-1) * log(tau)
                    - (1/tau) * (omegat_D * omega - (omegat_W * omega)));
  }

vector Lambda_bp(matrix G, vector gamma,int n){
      vector [n] lprob;
      lprob= -G*gamma;
      return lprob;
    }
vector log_lambda_bp(matrix g, vector gamma,int N,real zeta){
      vector[N] lprob1;
      lprob1 = log(g*gamma) - log(zeta);
            return lprob1;
    }
}

