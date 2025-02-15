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
