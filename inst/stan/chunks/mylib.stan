functions{

vector Lambda_plp1(vector max_stop, vector alpha, int n){
      vector [n] lprob;
      for(i in 1:n){
     lprob[i] =  -exp(lmultiply(alpha[2],max_stop[i]) - lmultiply(alpha[2],alpha[1]));
           }
      return lprob;
    }
vector log_lambda_plp1(vector time, int N, vector alpha){
      vector[N] lprob1;
      for(i in 1:N){
      lprob1[i] = log(alpha[2])-
      lmultiply(alpha[2], alpha[1])+
      lmultiply(alpha[2]-1,time[i]);
      // lprob1[i] = log(alpha[1])+log(alpha[2])+
      // lmultiply(alpha[2]-1,time[i]);
      }
      return lprob1;
    }

vector Lambda_plp2(vector max_stop, vector alpha, int n){
      vector [n] lprob;
      for(i in 1:n){
     // lprob[i] =  -exp(lmultiply(alpha[2],max_stop[i]) - lmultiply(alpha[2],alpha[1]));
      lprob[i]= -exp(log(alpha[1])+lmultiply(alpha[2],max_stop[i]));
      }
      return lprob;
    }
vector log_lambda_plp2(vector time, int N, vector alpha){
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

vector Lambda_plp3(vector max_stop, vector alpha, int n, real zeta){
      vector [n] lprob;
      vector [n] y1=max_stop/zeta;
    for(i in 1:n){
     lprob[i]= -exp(log(alpha[1])+lmultiply(alpha[2],y1[i]));
     // lprob[i] =  -exp(lmultiply(alpha[2],y1[i]) - lmultiply(alpha[2],alpha[1]));
           }
      return lprob;
    }
vector log_lambda_plp3(vector time, int N, vector alpha, real zeta){
      vector[N] lprob1;
      vector [N] y=time/zeta;
      for(i in 1:N){
      // lprob1[i] = log(alpha[2])-
      // lmultiply(alpha[2], alpha[1])+
      // lmultiply(alpha[2]-1,y[i]);
      lprob1[i] = log(alpha[1])+log(alpha[2])+
      lmultiply(alpha[2]-1,y[i]);
      }
      return (lprob1-log(zeta));
    }

vector Lambda_bp(matrix G, vector gamma,int n){
      vector [n] lprob;
      lprob= G*gamma;
      return lprob;
    }
vector log_lambda_bp(matrix g, vector gamma,int N,real zeta){
      vector[N] lprob1;
      lprob1 = log(g*gamma) - log(zeta);
            return lprob1;
    }
}
