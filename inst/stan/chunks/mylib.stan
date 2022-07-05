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
}
