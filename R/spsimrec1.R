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
spsimrec1 <- function(N,
                     x,
                     x1,
                     cov.rec,
                     cov.log = NULL,
                     beta.x.rec = NULL,
                     beta.x.log = NULL,
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
  dist.z <- tolower(dist.z)
  dist.z <- match.arg(dist.z)
  nr.cov.rec <- length(beta.x.rec)
  nr.cov.log <- length(beta.x.log)

  cov.rec<-x1 %>% select(cov.rec)
  cov.log<-x1 %>% select(cov.log)
  #print(cov.rec)
  #print(cov.log)

  ## gen_fu - Gera períodos de acompanhamento  ====
  # gen_fu<-function(N,cens.prob,fu.min,fu.max){
  #   fu <- rbinom(N, 1, cens.prob) # 1 = censored
  #   nr.cens <- sum(fu)
  #   if (nr.cens == 0) { # nobody censored
  #     fu <- runif(N, min = fu.min, max = fu.max)}
  #   else {
  #     index.cens <- which(fu == 1)
  #     fu[-index.cens] <- runif((N - nr.cens), min = fu.min, max = fu.max)
  #     fu[index.cens] <- runif(nr.cens, min = 0, max = fu.max)
  #   }
  #   return(fu)
  # }
  # set.seed(123)
  # fu<-gen_fu(N,cens.prob,fu.min,fu.max)
  # set.seed(NULL)
  #
  # ## gen_cov - Gera covariáveis  ====
  # gen_cov<-function(N,
  #                   dist.x,
  #                   par.x,
  #                   beta.x){
  #   nr.cov <- length(beta.x)
  #   x <- matrix(0, N, nr.cov)
  #   for (i in 1:nr.cov) {
  #     dist.x[i] <- match.arg(dist.x[i], choices = c("binomial", "normal"))
  #     if (dist.x[i] == "binomial") {
  #       x[, i] <- c(rbinom(N, 1, par.x[[i]]))
  #
  #     } else {
  #       mu.x <- par.x[[i]][1]
  #       sigma.x <- par.x[[i]][2]
  #       x[, i] <- c(rnorm(N, mean = mu.x, sd = sigma.x))
  #     }
  #   }
  #   return(x)
  # }
  #
  # set.seed(1)
  # if(nr.cov!=0){
  #   x<-gen_cov(N,
  #              dist.x,
  #              par.x,
  #              beta.x)
  #   x1<-as.data.frame(cbind(ID,x))
  #   colnames(x1)<-c("ID",paste("X",c(1:nr.cov),sep=""))
  # }else{x<-rep(0,N)
  # x1<-as.data.frame(cbind(ID,x))
  # colnames(x1)<-c("ID","X")}
  # set.seed(NULL)

  ## gen_rnd_ef - Gera efeitos aleatórios  ====
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

  if(random.ef!=0){
    #set.seed(321)
    rnd_ef<- gen_rnd_ef(N, ID, dist.z, tp_rnd_ef, par.z,mu.omega,sigma.omega)
    #set.seed(NULL)
  }else{
    rnd_ef<-rep(1,N)
    tp_rnd_ef<-0
  }

  rnd_ef1<-as.data.frame(cbind(ID,rnd_ef))
  #print(rnd_ef1)
  if(tp_rnd_ef==0){colnames(rnd_ef1)<-c("ID","z")}
  else{colnames(rnd_ef1)<-c("ID","w")}

  ## Define indivíduos recorrentes (INFLAÇÃO DE ZEROS)  ====

  gen_zi<-function(ID,N,pi,logist,cov.log,beta.x.log){
    if(logist==0){
  recurr <- rbinom(N, 1, pi)
  # recurr<-rbind(recurr,pi)
  pi<-rep(pi,N)
  return(list(recurr=recurr,pi=pi))
  }

  else {
   # pi<-1/(1+exp(-(1+as.matrix(cov.log) %*% beta.x.log)))
   pi<-1/(1+exp(-(as.matrix(cov.log) %*% beta.x.log)))
   recurr<-NULL
   for (i in 1:N){
     recurr[i] <- rbinom(1, 1,pi[i])
   }
   # recurr<-rbind(recurr,pi)
  }
    #print(recurr)
    #return(list(N=N,x=x,x1=x1,fu=fu,nr.cov=nr.cov,beta.x=beta.x,fu.max=fu.max,fu.min=fu.min))
    return(list(recurr=recurr,pi=pi))
    #return(recurr)
  }

  recurr<-gen_zi(ID,N,pi,logist,cov.log,beta.x.log)
  pi<-recurr$pi
  #print(pi)
  recurr<-recurr$recurr
  #print(recurr)
  recurr1<-as.data.frame((cbind(ID,recurr,pi)))

  colnames(recurr1)<-c("ID","recurr","pi")
 # print(recurr1)

#print(recurr1)

  ## gen_data - Gera tempos de ocorrência dos eventos  ====
  gen_data<-function(ID,
                     N,
                     dist.rec,
                     par.rec,
                     fu,
                     cov.rec,
                     beta.x.rec,
                     rnd_ef){

    if (dist.rec == "weibull") { # weibull
      alpha1 <- par.rec[1]
      alpha2 <- par.rec[2]
    }

    #print(cov.rec)
    #print(beta.x.rec)
    ## Cálculo de alpha1_este e exp_eta ====
    # Considera a forma utilizada para introdução de efeitos aleatórios
    if(nr.cov.rec==0){exp_eta=rep(1,N)}
    else if(tp_rnd_ef==0){#{Y_i(t) * \lambda_0(t)* Z_i *exp(\beta^t X_i)}
      exp_eta <- exp(as.matrix(cov.rec) %*% beta.x.rec) * rnd_ef
    }else{#{Y_i(t) * \lambda_0(t)*exp(\beta^t X_i+\omega_i)}
      exp_eta <- exp(as.matrix(cov.rec) %*% beta.x.rec + rnd_ef)
    }
    alpha1_eta <- alpha1*exp_eta


    ## Definição dos tempos de ocorrência dos primeiros eventos ====
    T<-NULL
    T1<-NULL
   # IND<-NULL
    for (i in 1:N) {
      t<-NULL
      U <- runif(1)
      if (dist.rec == "weibull") {
        t <- ((-1)*log(U)*(alpha1_eta[i])^(-1))^(1 / alpha2) # (veja artigo Generating survival times to simulate pag 1717 tabela II)
        #ind<-0
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
      select(-c(time,ngroup))
    return(tab)
    #set.seed(NULL)
  }
  tab <-gen_data(ID, N, dist.rec, par.rec, fu, cov.rec,beta.x.rec,rnd_ef)
  return(tab)
}
