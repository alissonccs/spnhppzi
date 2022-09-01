functions{
  //modelo sem fragilidade
vector Lambda_plp2(vector max_stop, vector alpha, int n){
      vector [n] lprob;
      for(i in 1:n){
     // lprob[i] =  -exp(lmultiply(alpha[2],max_stop[i]) - lmultiply(alpha[2],alpha[1]));
      lprob[i]= -exp(log(alpha[1])+lmultiply(alpha[2],max_stop[i]));
      }
      return lprob;
    }
vector lambda_plp2(vector time, int N, vector alpha){
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

real sparse_car_lpdf(vector phi, real tau, real alpha,
    int[,] W_sparse, vector D_sparse, vector lambda, int SP_N, int W_n) {
      row_vector[SP_N] phit_D; // phi' * D
      row_vector[SP_N] phit_W; // phi' * W
      vector[SP_N] ldet_terms;

      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, SP_N);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }

      for (i in 1:SP_N) ldet_terms[i] = log1m(alpha * lambda[i]);
      return 0.5 * (SP_N * log(tau)
                    + sum(ldet_terms)
                    - tau * (phit_D * phi - alpha * (phit_W * phi)));
  }

real sparse_iar_lpdf(vector phi, real tau,
    int[,] W_sparse, vector D_sparse, vector lambda, int SP_N, int W_n) {
      row_vector[SP_N] phit_D; // phi' * D
      row_vector[SP_N] phit_W; // phi' * W
      vector[SP_N] ldet_terms;

      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, SP_N);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }

      return 0.5 * ((SP_N-1) * log(tau)
                    - tau * (phit_D * phi - (phit_W * phi)));
  }
}
