#'@title spsimrec
#' @aliases spsimrec
#' @export
#' @description Recur:
# @param cens.prob  Probabiidade de apresentar censura antes do período máximo de acompanhamento.
# @param dist.x     Distribuição das Covariáveis. Binomial ou normal.
# @param par.x      Parâmetros das distribuições das covariáveis.
#' @param N          Número de indivíduos
#' @param spatial    Indicadora do uso de modelo espacial, (default =0 , não espacial)
#' @param sp_model   Tipo de modelo espacial
#' @param SP_N       Número de unidades de área
#' @param nb_mat     Matriz de vizinhança
#' @param sp_tau     Valor do parâmetro de precisão no modelo ICAR (\tau)
#' @param sp_alpha   Valor do associação espacial no modelo CAR.
#' @param x          Matriz de covariávies
#' @param x1         Matriz ID e covariávies
#' @param cov_rec    Lista de covariáveis relacionadas à função de intensidade
#' @param cov_log    Lista de covariáveis relacionadas à regressão logística
#' @param beta_x_rec Coeficientes das covariáveis relacionadas à função de intensidade
#' @param beta_x_log Coeficientes das covariáveis relacionadas à regressão logística
#' @param xi         Parâmetro que relaciona a função de intensidade e a logística
#' @param fu_min     Tempo mínimo de acompanhamento
#' @param fu_max     Tempo máximo de acompanhamento
#' @param random_ef  Indicadora de uso de efeitos aleatórios na geração dos dados (random_ef==0, modelo sem efeitos).
#' @param dist_z     Distribuição do efeito aleatório não espacial (Gamma ou lognormal).
#' @param tp_rnd_ef  Indicadora de efeito aleatório multiplicativo
#' @param par_z      Parâmetros de variância da distribuição do efeito aleatório.
#' @param dist_int_func   Forma da função de intensidade.  "Weibull" (Lei de potência)
#' @param par_int_func    Parâmetros da função de intensidade. Escala e forma (Lei de potência)
#' @param pi         Inflação de zeros (proporção de indivíduos sem recorrências)
#' @param logist     Indicadora do uso de covariáveis na logística

# SIMRECEV - SIMULAÇÃO DE EVENTOS RECORRENTES ====
spsimrec7 <-  function(N,
                      spatial = 0,
                      sp_model = c("car","sparse","icar"),
                      list_area=NULL,
                      SP_N=NULL,
                      nb_mat=NULL,
                      sp_tau=NULL,
                      sp_alpha=NULL,
                      # beta.x = NULL,
                      # x,
                      x1,
                      cov_rec = NULL,
                      cov_log = NULL,
                      beta_x_rec = NULL,
                      beta_x_log = NULL,
                      xi=1,
                      fu,
                      fu_max,
                      fu_min,
                      random_ef=0,
                      dist_z = c("gamma","lognormal"),
                      tp_rnd_ef=0,
                      par_z=0,
                      dist_int_func,
                      par_int_func=0,
                      pi = 0,
                      # dfree = 0,
                      logist = 0,
                      mu_omega=0,
                      sigma_omega=0){
  ID <- c(1:N)
  dist_z <- tolower(dist_z)
  dist_z <- match.arg(dist_z)
  nr.cov_rec <- length(beta_x_rec)
  nr.cov_log <- length(beta_x_log)
  cov_rec<-x1 %>% dplyr::select(all_of(cov_rec))
  cov_log<-x1 %>% dplyr::select(all_of(cov_log))

  ## Designa as unidades de areas para cada individuo a partir de uma lista de areas fornecida como entrada ====
  if(spatial==1){
  SP_ID<- sample(list_area, N, replace=TRUE, prob=rep(1/SP_N,SP_N) )
  ID_SP_ID<-as.data.frame(cbind(ID,SP_ID))
  }

  #nr.cov <- length(beta.x)
  sp_model<-tolower(sp_model)
  sp_model<-match.arg(sp_model)

  sp_model<-switch(sp_model,
                   "car" = 1,
                   "sparse"=2,
                   "icar" = 3)

  # EFEITOS ALEATÓRIOS  ====
  ### Funções ====
  #### Efeito não espacial ====
  gen_rnd_ef<-function(N, ID, dist_z, tp_rnd_ef, par_z,mu_omega,sigma_omega){
    if (tp_rnd_ef==0){ #Entra com parâmetros para Z. {Y_i(t) * \lambda_0(t)* Z_i *exp(\beta^t X_i)}
      if(par_z==0){# se par_z=0 então a fragilidade=1 para todos
        z <- rep(1, N)}
      else{
        dist_z <- match.arg(dist_z, choices = c("gamma", "lognormal"))
        if (dist_z == "gamma") { # gamma-frailty
          aGamma <- 1 / par_z
          rnd_ef <- rgamma(N, shape = aGamma, scale = 1 / aGamma)}
        else { # lognormal  -- E(Z)=1
          mu <- log(1 / sqrt(par_z + 1))
          sigma <- sqrt(log(par_z + 1))
          rnd_ef <- exp(rnorm(N, mean = mu, sd = sigma))
        }
      }
    }
    else{ #Entra com parâmetros para \omega. {Y_i(t) * \lambda_0(t)*exp(\beta^t X_i+\omega_i)}
      rnd_ef <- rnorm(N, mean = mu_omega, sd = sigma_omega)
    }
    return(rnd_ef)
  }

  #### Efeito espacial ====
  ##### CAR ====
  car_sp_rnd_ef<-function(SP_N, sp_tau, sp_alpha, nb_mat){
      print("CAR")
      D<-diag(rowSums(nb_mat))
      prec_mat<-sp_tau*(D-sp_alpha*nb_mat)
      #rnd_ef <- mvrnorm(n = 1, mu =rep(0,SP_N), Sigma =(SIGMA))
     rnd_ef <- mvrnorm(n = 1, mu =rep(0,SP_N), Sigma = solve(prec_mat))
      print(paste("rnd_ef: ", rnd_ef))
      print(paste("rnd_ef mean: ", mean(rnd_ef)))
    return(rnd_ef)
  }

  ##### CAR CHOLESKY ====
  CAR.simWmat <- function(sp_tau, sp_alpha, nb_mat){
    print("CAR SPARSE")
    D<-diag(rowSums(nb_mat))
    B<-solve(D)%*%nb_mat
    n <- nrow(B)
    I <- diag.spam(1, n)
    Q <- as.spam(sp_tau*D%*% (I - sp_alpha*B))
    ## Simulate from N(0, Q)
    cholR <- chol.spam(Q, pivot = "MMD", memory = list(nnzcolindices = 6.25 * n)) ## upper triagle
    rnd_ef <- backsolve(cholR, rnorm(n))
    ## return results
    return(rnd_ef)
  }

  ##### ICAR ====
   icar_sp_rnd_ef <- function(nb_mat,sig=1){
    print("ICAR")
    num <- rowSums(nb_mat)
    n <- ncol(nb_mat)

    Q <- -nb_mat
    diag(Q) <- num

    Q_aux=eigen(Q)$vectors[,order(eigen(Q)$values)]

    D_aux=sort(eigen(Q)$values)

    rnd_ef <- rnorm(n-1,0,sqrt(sig*(1/D_aux[-1])))
    rnd_ef <- Q_aux%*%c(0,rnd_ef)
    return(as.vector(rnd_ef))
  }


### Executa funções para gerar efeitos aleatórios ====
   if(random_ef==0){  #sem efeitos
     rnd_ef<-rep(1,N)
     tp_rnd_ef<-0
     rnd_ef1<-as.data.frame(cbind(ID,rnd_ef))
     colnames(rnd_ef1)<-c("ID","rnd_ef")
   }else if(spatial==0){ #efeitos não espaciais
     rnd_ef<- gen_rnd_ef(N, ID, dist_z, tp_rnd_ef, par_z,mu_omega,sigma_omega)
     rnd_ef1<-as.data.frame(cbind(ID,rnd_ef))

     if(tp_rnd_ef==0){colnames(rnd_ef1)<-c("ID","rnd_ef")}
     else{colnames(rnd_ef1)<-c("ID","rnd_ef")}
   }else{  #efeitos espaciais
     tp_rnd_ef<-1
     if (sp_model==1){ #CAR
       rnd_ef<-car_sp_rnd_ef(SP_N=SP_N, sp_tau=sp_tau,sp_alpha=sp_alpha, nb_mat=nb_mat)
     }
     if (sp_model==2){ #CAR CHOLESKY
       rnd_ef<-CAR.simWmat(sp_tau=sp_tau,sp_alpha = sp_alpha,nb_mat=nb_mat)
     }
     if (sp_model==3){ #ICAR
       rnd_ef<-icar_sp_rnd_ef(nb_mat,sig=1/sp_tau)
     }

     rnd_ef<-as.data.frame(cbind(as.numeric(row.names(nb_mat)),rnd_ef))
     colnames(rnd_ef)<-c("SP_ID","rnd_ef")
     rnd_ef1<-ID_SP_ID %>%
       left_join(rnd_ef,by="SP_ID")
     colnames(rnd_ef1)<-c("ID","SP_ID","rnd_ef")
     # rnd_ef<-rnd_ef$rnd_ef
   }

   # INFLAÇÃO DE ZEROS - Define indivíduos recorrentes  ====
   ### Função ====
   gen_zi<-function(ID,N,pi,logist,cov_log,beta_x_log,xi1=1,rnd_ef){
     if(logist==0){
       recurr <- rbinom(N, 1, pi)
       pi<-rep(pi,N)
       # return(list(recurr=recurr,pi=pi))
       #recurr<-rbind(recurr,pi)
       # pi<-rep(pi,N)
       # return(recurr)
     }
     else {
       pi<-1/(1+exp(-(1+as.matrix(cov_log) %*% beta_x_log)))
       if(is.null(xi) == TRUE){
         pi<-1/(1+exp(-(as.matrix(cov_log) %*% beta_x_log +rnd_ef)))
       }else{
         pi<-1/(1+exp(-(as.matrix(cov_log) %*% beta_x_log + xi1*rnd_ef)))
       }
       recurr<-NULL
       for (i in 1:N){
         recurr[i] <- rbinom(1, 1,pi[i])
       }
       # recurr<-rbind(recurr,pi)
     }
     #print(recurr)
     #return(list(N=N,x=x,x1=x1,fu=fu,nr.cov=nr.cov,beta.x=beta.x,fu_max=fu_max,fu_min=fu_min))
     return(list(recurr=recurr,pi=pi))
     #return(recurr)
   }
 ### Executa função ====
   recurr<-gen_zi(ID,N,pi,logist,cov_log,beta_x_log,xi,rnd_ef1$rnd_ef)
   pi<-recurr$pi
   #print(pi)
   recurr<-recurr$recurr
   #print(recurr)
   recurr1<-as.data.frame((cbind(ID,recurr,pi)))

   colnames(recurr1)<-c("ID","recurr","pi")

   input_gen_data<-cov_rec %>%
     left_join(rnd_ef1,by="ID") %>%
     left_join(recurr1,by="ID") %>%
     left_join(fu,by="ID")


  ## TEMPO DE OCORRÊNCIA DOS EVENTOS ====
   gen_data<-function(ID,
                      N,
                      dist_int_func,
                      par_int_func,
                      fu,
                      x,
                      rnd_ef,
                      recurr){

    if (dist_int_func == "weibull") { # weibull
      alpha1 <- par_int_func[1]
      alpha2 <- par_int_func[2]
    }

    ## Cálculo de alpha1_este e exp_eta ====
    # Considera a forma utilizada para introdução de efeitos aleatórios
    if(nr.cov_rec==0){exp_eta=rep(1,N)}
    else if(tp_rnd_ef==0){#{Y_i(t) * \lambda_0(t)* Z_i *exp(\beta^t X_i)}
      exp_eta <- exp(x %*% beta_x_rec) * rnd_ef
    }else{#{Y_i(t) * \lambda_0(t)*exp(\beta^t X_i+\omega_i)}
      exp_eta <- exp(x %*% beta_x_rec + rnd_ef)
    }
    alpha1_eta <- alpha1*exp_eta
   # print(paste("alpha1_eta: ", dim(alpha1_eta)))
   # print(paste("exp_eta: ", dim(exp_eta)))
   # print(paste("N: ", N))
    #print(alpha2)

    ## Definição dos tempos de ocorrência dos primeiros eventos ====

    T<-NULL
    T1<-NULL
   # print(paste("T: ", dim(T)))
   # print(paste("T1: ", dim(T1)))
   # IND<-NULL
    for (i in 1:N) {
      t<-NULL
      #print(t)
      U <- runif(1)
      if (dist_int_func == "weibull") {
        t <- ((-1)*log(U)*(alpha1_eta[i])^(-1))^(1 / alpha2) # (veja artigo Generating survival times to simulate pag 1717 tabela II)
        #ind<-0
        # print(i)
        # print(t)
         # print(fu[i])
        if (t>fu[i]){
          t<-fu[i]
         # ind<-1
          }
      }
      T1 <- cbind(ID[i],t)

      ## Definição dos tempos de ocorrência dos eventos subsequentes ====
      if (recurr[i]==0 & t<fu[i]){
        # print(ID[i])
        # print(recurr[i])
        while (t < fu[i]) {
          U <- runif(1)
          t1 <- t
          if (dist_int_func == "weibull") { # weibull
            t <- ((-1)*log(U)*(alpha1_eta[i])^(-1) + (t1)^(alpha2))^(1 / alpha2)
          }
          #print(t)
          if (t >= fu[i]) break
          T1 <- rbind(T1,c(ID[i],t))
         # print(T1)
        }
      }
      T<-as.data.frame(rbind(T,T1))
      # print(T)
    }
    colnames(T)<-c("ID","time")
    # print(paste("T: ", dim(T)))

    ## Consolida tabela contendo os dados de saída ====
    tab <-T %>%
      group_by(ID)%>%
      mutate(#individuo = group_indices(),
        ngroup=n(),
        rep=row_number(),
        expand=case_when((ngroup==rep & !(rep==1&(time==0 | time==fu_max)))~2,TRUE~1),
        expand1=expand)%>%
      expandRows("expand") %>%
      mutate(ngroup1=n(),
             IndRec=case_when(ngroup1>2~1, TRUE~0),
             rep1=row_number(),
             begin=case_when(rep1==1~0,TRUE~lag(time)),
             end=case_when(ngroup1==rep1~fu_max, TRUE~time),
             status=case_when(end==fu_max~0,TRUE~1)) %>%
      ungroup() %>%
      left_join(x1,by="ID") %>%
      left_join(rnd_ef1,by="ID") %>%
      left_join(recurr1,by="ID") %>%
      dplyr::select(-c(time,ngroup))
    return(tab)
    #set.seed(NULL)
   }
   tab <-gen_data(ID=input_gen_data$ID,N=N, dist_int_func=dist_int_func, par_int_func=par_int_func, fu=input_gen_data$fu, x=as.matrix(input_gen_data[,2:(1+nr.cov_rec)]),rnd_ef=input_gen_data$rnd_ef,input_gen_data$recurr)
  #tab <-gen_data(ID, N, dist_int_func, par_int_func, fu, x,rnd_ef)
  return(tab)
}
