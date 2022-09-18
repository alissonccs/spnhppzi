#'@title spsimrec
#' @aliases spsimrec
#' @export
#' @description Recur:
# @param cens.prob  Probabiidade de apresentar censura antes do período máximo de acompanhamento.
# @param dist.x     Distribuição das Covariáveis. Binomial ou normal.
# @param par.x      Parâmetros das distribuições das covariáveis.
#' @param beta.x     Coeficiente de regressão das covariáveis.
#' @param dist.z     Distribuição do efeito aleatório. Gamma ou lognormal.
#' @param par.z      Parâmetros da distribuição do efeito aleatório.
#' @param dist.rec   Forma da função de intensidade.  "Weibull" (Lei de potência)
#' @param par.rec    Parâmetros da função de intensidade. Escala e forma

# SIMRECEV - SIMULAÇÃO DE EVENTOS RECORRENTES ====
spsimrec3 <-  function(N,
                      nr.cov,
                      spatial,
                      sp_model = c("car","sparse","icar"),
                      list_area,
                      SP_N,
                      nb_mat,
                      sp_tau,
                      sp_alpha,
                      beta.x = NULL,
                      x,
                      x1,
                      fu,
                      fu.max,
                      fu.min,
                      dist.z = c("gamma","lognormal"),
                      random.ef=0,
                      tp_rnd_ef=0,
                      par.z=0,
                      dist.rec,
                      par.rec=0,
                      pi = 0,
                      dfree = 0,
                      logist = 0,
                      mu.omega=0,
                      sigma.omega=0){
  ID <- c(1:N)
  #SP_ID<-c(1:SP_N)
  # list_RMBH<-c(3165537,3141108,3126000,3140704,3137601,3136652,3132206,3130101,
  #              3129806,3157807,3155306,3154804,3154606,3153905,3124104,3171204,
  #              3168309,3162955,3156700,3133709,3134608,3136603,3140159,3144805,
  #              3106705,3110004,3149309,3105004,3106200,3109006,3112505,3117876,
  #              3118601, 3162922
  # )
  # N<-1000
  # SP_N<-34

  if(spatial==1){
  SP_ID<- sample( list_area, N, replace=TRUE, prob=rep(1/SP_N,SP_N) )
  }
  ID_SP_ID<-as.data.frame(cbind(ID,SP_ID))
 # ID_SP_ID<-cbind(ID,SP_ID)
  dist.z <- tolower(dist.z)
  dist.z <- match.arg(dist.z)
  #nr.cov <- length(beta.x)
  sp_model<-tolower(sp_model)
  sp_model<-match.arg(sp_model)

  sp_model<-switch(sp_model,
                   "car" = 1,
                   "sparse"=2,
                   "icar" = 3)

  ## gen_rnd_ef - Gera efeitos aleatórios  ====

  # Efeito não espacial
  gen_rnd_ef<-function(N, ID, dist.z, tp_rnd_ef, par.z,mu.omega,sigma.omega){
    if (tp_rnd_ef==0){ #Entra com parâmetros para Z. {Y_i(t) * \lambda_0(t)* Z_i *exp(\beta^t X_i)}
      if(par.z==0){# se par.z=0 então a fragilidade=1 para todos
        z <- rep(1, N)}
      else{
        dist.z <- match.arg(dist.z, choices = c("gamma", "lognormal"))
        if (dist.z == "gamma") { # gamma-frailty
          aGamma <- 1 / par.z
          rnd_ef <- rgamma(N, shape = aGamma, scale = 1 / aGamma)}
        else { # lognormal  -- E(Z)=1
          mu <- log(1 / sqrt(par.z + 1))
          sigma <- sqrt(log(par.z + 1))
          rnd_ef <- exp(rnorm(N, mean = mu, sd = sigma))
        }
      }
    }
    else{ #Entra com parâmetros para \omega. {Y_i(t) * \lambda_0(t)*exp(\beta^t X_i+\omega_i)}
      rnd_ef <- rnorm(N, mean = mu.omega, sd = sigma.omega)
    }
    return(rnd_ef)
  }

  # Efeito espacial
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



   if(random.ef==0){
     rnd_ef<-rep(1,N)
     tp_rnd_ef<-0
     colnames(rnd_ef1)<-c("ID","rnd_ef")
   }else if(spatial==0){
     rnd_ef<- gen_rnd_ef(N, ID, dist.z, tp_rnd_ef, par.z,mu.omega,sigma.omega)
     rnd_ef1<-as.data.frame(cbind(ID,rnd_ef))

     if(tp_rnd_ef==0){colnames(rnd_ef1)<-c("ID","rnd_ef")}
     else{colnames(rnd_ef1)<-c("ID","rnd_ef")}
   }else{
     tp_rnd_ef<-1
     if (sp_model==1){
       rnd_ef<-car_sp_rnd_ef(SP_N=SP_N, sp_tau=sp_tau,sp_alpha=sp_alpha, nb_mat=nb_mat)
     }
     if (sp_model==2){
       rnd_ef<-CAR.simWmat(sp_tau=sp_tau,sp_alpha = sp_alpha,nb_mat=nb_mat)
     }
     if (sp_model==3){
       rnd_ef<-icar_sp_rnd_ef(nb_mat,sig=1/sp_tau)
     }

     rnd_ef<-as.data.frame(cbind(as.numeric(row.names(nb_mat)),rnd_ef))
     colnames(rnd_ef)<-c("SP_ID","rnd_ef")
     rnd_ef1<-ID_SP_ID %>%
       left_join(rnd_ef,by="SP_ID")
     colnames(rnd_ef1)<-c("ID","SP_ID","rnd_ef")
     # rnd_ef<-rnd_ef$rnd_ef
   }


  ## Define indivíduos recorrentes (INFLAÇÃO DE ZEROS)  ====

   gen_zi<-function(ID,N,pi){
     recurr <- t(rbinom(N, 1, pi))
   }
   #set.seed(234)
   # recurr<-gen_zi(ID,N,pi)
   # pi<-0
   recurr<-gen_zi(ID,N,pi)
   # set.seed(NULL)

   recurr1<-as.data.frame(t(rbind(ID,recurr)))
   colnames(recurr1)<-c("ID","recurr")

   # if(logist==1){
   #  pi<-1/(1+exp(-(1+x %*% beta.x)))
   # }


   input_gen_data<-x1 %>%
     left_join(rnd_ef1,by="ID") %>%
     left_join(recurr1,by="ID") %>%
     left_join(fu,by="ID")


  ## gen_data - Gera tempos de ocorrência dos eventos  ====
   gen_data<-function(ID,
                      N,
                      dist.rec,
                      par.rec,
                      fu,
                      x,
                      rnd_ef,
                      recurr){

    if (dist.rec == "weibull") { # weibull
      alpha1 <- par.rec[1]
      alpha2 <- par.rec[2]
    }

    ## Cálculo de alpha1_este e exp_eta ====
    # Considera a forma utilizada para introdução de efeitos aleatórios
    if(nr.cov==0){exp_eta=rep(1,N)}
    else if(tp_rnd_ef==0){#{Y_i(t) * \lambda_0(t)* Z_i *exp(\beta^t X_i)}
      exp_eta <- exp(x %*% beta.x) * rnd_ef
    }else{#{Y_i(t) * \lambda_0(t)*exp(\beta^t X_i+\omega_i)}
      exp_eta <- exp(x %*% beta.x + rnd_ef)
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
      if (dist.rec == "weibull") {
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
          if (dist.rec == "weibull") { # weibull
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
    print(paste("T: ", dim(T)))

    ## Consolida tabela contendo os dados de saída ====
    tab <-T %>%
      group_by(ID)%>%
      mutate(#individuo = group_indices(),
        ngroup=n(),
        rep=row_number(),
        expand=case_when((ngroup==rep & !(rep==1&(time==0 | time==fu.max)))~2,TRUE~1),
        expand1=expand)%>%
      expandRows("expand") %>%
      mutate(ngroup1=n(),
             IndRec=case_when(ngroup1>2~1, TRUE~0),
             rep1=row_number(),
             begin=case_when(rep1==1~0,TRUE~lag(time)),
             end=case_when(ngroup1==rep1~fu.max, TRUE~time),
             status=case_when(end==fu.max~0,TRUE~1)) %>%
      ungroup() %>%
      left_join(x1,by="ID") %>%
      left_join(rnd_ef1,by="ID") %>%
      left_join(recurr1,by="ID") %>%
      dplyr::select(-c(time,ngroup))
    return(tab)
    #set.seed(NULL)
   }
   tab <-gen_data(ID=input_gen_data$ID,N=N, dist.rec=dist.rec, par.rec=par.rec, fu=input_gen_data$fu, x=as.matrix(input_gen_data[,2:(1+nr.cov)]),rnd_ef=input_gen_data$rnd_ef,input_gen_data$recurr)
  #tab <-gen_data(ID, N, dist.rec, par.rec, fu, x,rnd_ef)

  return(tab)
}
