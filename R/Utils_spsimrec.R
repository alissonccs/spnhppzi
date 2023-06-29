  ### Funções ====
  # EFEITOS ALEATÓRIOS  ----
  #### Efeito não espacial ====
  #' @title gen_rnd_ef_
  #' @aliases gen_rnd_ef
  #' @export
  #' @description Gera efeitos aleatórios não espacial.
  #' @param N          Número de indivíduos.
  #' @param ID         Identificador de indivíduos.
  #' @param random_ef  Indicadora de uso de efeitos aleatórios na geração dos dados (random_ef==0, modelo sem efeitos).
  #' @param dist_z     Distribuição do efeito aleatório não espacial (Gamma ou lognormal).
  #' @param tp_rnd_ef  Indicadora de efeito aleatório multiplicativo.
  #' @param par_z      Parâmetros de variância da distribuição do efeito aleatório.
  #' @param mu_omega   Parâmetro média da distribuição dos efeitos aleatórios para a forma {Y_i(t)*lambda_0(t)*exp(beta^t X_i+omega_i)}
  #' @param sigma_omega   Parâmetro desvio padrão da distribuição dos efeitos aleatórios para a forma {Y_i(t)*lambda_0(t)*exp(beta^t X_i+omega_i)}

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
  #' @title car_sp_rnd_ef
  #' @aliases car_sp_rnd_ef
  #' @export
  #' @description      Gera efeitos aleatórios espaciais considerando a estrutura do modelo CAR.
  #' @param SP_N       Número de unidades de área.
  #' @param nb_mat     Matriz de vizinhança.
  #' @param sp_tau     Valor do parâmetro de precisão no modelo ICAR (\tau).
  #' @param sp_alpha   Valor do associação espacial no modelo CAR.

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
  #' @title CAR.simWmat
  #' @aliases CAR.simWmat
  #' @export
  #' @description      Gera efeitos aleatórios espaciais utilizando a decomposição de Cholesky.
  #' @param nb_mat     Matriz de vizinhança.
  #' @param sp_tau     Valor do parâmetro de precisão no modelo ICAR (\tau).
  #' @param sp_alpha   Valor do associação espacial no modelo CAR.

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
  #' @title icar_sp_rnd_ef
  #' @aliases icar_sp_rnd_ef
  #' @export
  #' @description      Gera efeitos aleatórios espaciais utilizando a estrutura do modelo ICAR.
  #' @param nb_mat     Matriz de vizinhança.
  #' @param sig        Parâmetro desvio padrão dos efeitos aleatórios.

  icar_sp_rnd_ef <- function(nb_mat,sig=1){
    print("ICAR")
    # print(head(nb_mat))
    num <- rowSums(nb_mat)
    n <- ncol(nb_mat)

    Q <- -nb_mat
    diag(Q) <- num

    Q_aux=eigen(Q)$vectors[,order(eigen(Q)$values)]

    D_aux=sort(eigen(Q)$values)

    rnd_ef <- rnorm(n-1,0,sqrt(sig*(1/D_aux[-1])))
    rnd_ef <- Q_aux%*%c(0,rnd_ef)
    # print((rnd_ef))
    return(as.vector(rnd_ef))
  }





  # INFLAÇÃO DE ZEROS - Define indivíduos recorrentes  ====
  #' @title spsimrec
  #' @aliases spsimrec
  #' @export
  #' @description Recur:
  #' @param N          Número de indivíduos.
  #' @param cov_log    Lista de covariáveis relacionadas à regressão logística.
  #' @param beta_x_log Coeficientes das covariáveis relacionadas à regressão logística.
  #' @param xi         Parâmetro que relaciona a função de intensidade e a logística.
  #' @param rnd        Vetor de efeitos aleatórios.
  #' @param pi_zi         Inflação de zeros (proporção de indivíduos sem recorrências).

  gen_zi<-function(ID,N,pi_zi,logist,cov_log,beta_x_log,xi1=1,rnd_ef){
    if(logist==0){ #Sem covariáveis na regressão logística - propabilidade de recorrência comum a todos os indivíduos
      recurr <- rbinom(N, 1, pi_zi)
      pi_zi<-rep(pi_zi,N)
      }
    else { # Com covariáveis na regressão logística - cada indivíduo possui sua probabilidade de recorrência definida em função das covariáveis e efeitos aleatórios.
      pi_zi<-1/(1+exp(-(1+as.matrix(cov_log) %*% beta_x_log)))
      # if(is.null(xi1) == TRUE){
      #   pi_zi<-1/(1+exp(-(as.matrix(cov_log) %*% beta_x_log +rnd_ef)))
      # }else{
      #   pi_zi<-1/(1+exp(-(as.matrix(cov_log) %*% beta_x_log + xi1*rnd_ef)))
      # }
      recurr<-NULL
      for (i in 1:N){
        recurr[i] <- rbinom(1, 1,pi_zi[i])
      }
    }
    return(list(recurr=recurr,pi_zi=pi_zi))
  }

  ## TEMPO DE OCORRÊNCIA DOS EVENTOS ====
  # @title gen_data
  # @aliases gen_data
  # @export
  # @description           Gera os tempos de recorrências dos eventos.
  # @param ID              Identificador dos indivíduos.
  # @param N               Número de indivíduos.
  # @param dist_int_func   Forma da função de intensidade.  "Weibull" (Lei de potência)
  # @param par_int_func    Parâmetros da função de intensidade. Escala e forma (Lei de potência)
  # @param rnd_ef          Vetor de efeitos aleatórios.
  # @param x               Matriz de covariávies
  # @param recurr          Indicador de recorrência.
  #
  # gen_data<-function(ID,
  #                    N,
  #                    dist_int_func,
  #                    par_int_func,
  #                    fu,
  #                    fu_max,
  #                    x,
  #                    # x1,
  #                    beta_x_rec,
  #                    tp_rnd_ef,
  #                    rnd_ef,
  #                    # rnd_ef1,
  #                    recurr,
  #                    # recurr1,
  #                    nr.cov_rec){
  #
  #   if (dist_int_func == "weibull") { # weibull
  #     alpha1 <- par_int_func[1]
  #     alpha2 <- par_int_func[2]
  #   }
  #
  #   ## Cálculo de alpha1_este e exp_eta ====
  #   # Considera a forma utilizada para introdução de efeitos aleatórios
  #   if(nr.cov_rec==0){exp_eta=rep(1,N)}
  #   else if(tp_rnd_ef==0){#{Y_i(t) * \lambda_0(t)* Z_i *exp(\beta^t X_i)}
  #     exp_eta <- exp(x %*% beta_x_rec) * rnd_ef
  #   }else{#{Y_i(t) * \lambda_0(t)*exp(\beta^t X_i+\omega_i)}
  #     exp_eta <- exp(x %*% beta_x_rec + rnd_ef)
  #   }
  #   alpha1_eta <- alpha1*exp_eta
  #   # print(paste("alpha1_eta: ", dim(alpha1_eta)))
  #   # print(paste("exp_eta: ", dim(exp_eta)))
  #   # print(paste("N: ", N))
  #   #print(alpha2)
  #
  #   ## Definição dos tempos de ocorrência dos primeiros eventos ====
  #
  #   T<-NULL
  #   T1<-NULL
  #   # print(paste("T: ", dim(T)))
  #   # print(paste("T1: ", dim(T1)))
  #   # IND<-NULL
  #   for (i in 1:N) {
  #     t<-NULL
  #     #print(t)
  #     U <- runif(1)
  #     if (dist_int_func == "weibull") {
  #       t <- ((-1)*log(U)*(alpha1_eta[i])^(-1))^(1 / alpha2) # (veja artigo Generating survival times to simulate pag 1717 tabela II)
  #       #ind<-0
  #       # print(i)
  #       # print(t)
  #       # print(fu[i])
  #       if (t>fu[i]){
  #         t<-fu[i]
  #         # ind<-1
  #       }
  #     }
  #     T1 <- cbind(ID[i],t)
  #
  #     ## Definição dos tempos de ocorrência dos eventos subsequentes ====
  #     if (recurr[i]==0 & t<fu[i]){
  #       # print(ID[i])
  #       # print(recurr[i])
  #       while (t < fu[i]) {
  #         U <- runif(1)
  #         t1 <- t
  #         if (dist_int_func == "weibull") { # weibull
  #           t <- ((-1)*log(U)*(alpha1_eta[i])^(-1) + (t1)^(alpha2))^(1 / alpha2)
  #         }
  #         #print(t)
  #         if (t >= fu[i]) break
  #         T1 <- rbind(T1,c(ID[i],t))
  #         # print(T1)
  #       }
  #     }
  #     T<-as.data.frame(rbind(T,T1))
  #     # print(T)
  #   }
  #   colnames(T)<-c("ID","time")
  #   # print(paste("T: ", dim(T)))
  #
  #   ## Consolida tabela contendo os dados de saída ====
  #   tab <-T %>%
  #     group_by(ID)%>%
  #     mutate(#individuo = group_indices(),
  #       ngroup=n(),
  #       rep=row_number(),
  #       expand=case_when((ngroup==rep & !(rep==1&(time==0 | time==fu_max)))~2,TRUE~1),
  #       expand1=expand)%>%
  #     expandRows("expand") %>%
  #     mutate(ngroup1=n(),
  #            IndRec=case_when(ngroup1>2~1, TRUE~0),
  #            rep1=row_number(),
  #            begin=case_when(rep1==1~0,TRUE~lag(time)),
  #            end=case_when(ngroup1==rep1~fu_max, TRUE~time),
  #            status=case_when(end==fu_max~0,TRUE~1)) %>%
  #     ungroup() %>%
  #     # left_join(x1,by="ID") %>%
  #     # left_join(rnd_ef1,by="ID") %>%
  #     # left_join(recurr1,by="ID") %>%
  #     dplyr::select(-c(time,ngroup))
  #   return(tab)
  #   #set.seed(NULL)
  # }



  ## TEMPO DE OCORRÊNCIA DOS EVENTOS ====
  #' @title gen_data_plp1
  #' @aliases gen_data_plp1
  #' @export
  #' @description           Gera os tempos de recorrências dos eventos.
  #' @param ID              Identificador dos indivíduos.
  #' @param N               Número de indivíduos.
  #' @param dist_int_func   Forma da função de intensidade.  "Weibull" (Lei de potência)
  #' @param par_int_func    Parâmetros da função de intensidade. Escala e forma (Lei de potência)
  #' @param rnd_ef          Vetor de efeitos aleatórios.
  #' @param x               Matriz de covariávies
  #' @param recurr          Indicador de recorrência.
  #'
  gen_data<-function(ID,
                     N,
                     dist_int_func,
                     par_int_func,
                     fu,
                     fu_max,
                     x,
                     # x1,
                     beta_x_rec,
                     tp_rnd_ef,
                     rnd_ef,
                     # rnd_ef1,
                     recurr,
                     # recurr1,
                     nr.cov_rec,
                     baseline){

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
    if (baseline==2){
    alpha1_eta <- alpha1*exp_eta
    }
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
        if(baseline==1){
        t <- alpha1*((-1)*log(U)*(exp_eta[i])^(-1))^(1 / alpha2) # (veja artigo Generating survival times to simulate pag 1717 tabela II)
        }
        else {#baseline==2
        t <- ((-1)*log(U)*(alpha1_eta[i])^(-1))^(1 / alpha2)
        }
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
            if(baseline==1){
            t <- ((-1)*log(U)*((exp_eta[i])^(-1))*alpha1^alpha2 + (t1)^(alpha2))^(1 / alpha2)
            }
            else{#baseline==2
            t <- ((-1)*log(U)*(alpha1_eta[i])^(-1) + (t1)^(alpha2))^(1 / alpha2)
            }
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
      # left_join(x1,by="ID") %>%
      # left_join(rnd_ef1,by="ID") %>%
      # left_join(recurr1,by="ID") %>%
      dplyr::select(-c(time,ngroup))
    return(tab)
    #set.seed(NULL)
  }


  ## TEMPO DE OCORRÊNCIA DOS PRIMIEROS EVENTOS - FUNC INTENSIDADE POLINOMIAL ====
  #' @title get_nhpp_realization_0
  #' @aliases get_nhpp_realization_0
  #' @export
  #' @description           Gera os tempos de ocorrência dos primeiros eventos.
  #' @param lambda          Funćão de intensidade polinomial.

  get_nhpp_realization_0 <- function(lambda,fu_max, x, beta_x_rec,rnd_ef){
    t_max <- fu_max
    t <- 0
    s <- 0
    Lambda <- function(tupper) integrate(f = lambda,x=x,beta_x_rec=beta_x_rec, rnd_ef, lower = 0, upper = tupper)$value
    Lambda_inv <- function(s){
      v <- seq(0,t_max+3, length.out = 250)
      min(v[Vectorize(Lambda)(v)>=s])
    }

    X <- numeric(0)
    # while(t <= t_max){
    u <- runif(1)
    s <- s -log(u)
    t <- Lambda_inv(s)
    X <- c( X, t)
    # }

    return(X)
  }



  ## TEMPO DE OCORRÊNCIA DAS RECORRÊNCIAS - FUNC INTENSIDADE POLINOMIAL ====
  #' @title get_nhpp_realization
  #' @aliases get_nhpp_realization
  #' @export
  #' @description           Gera os tempos de ocorrência das reincidências.
  #' @param lambda          Funćão de intensidade polinomial.


  get_nhpp_realization <- function(lambda,fu_max, x, beta_x_rec,rnd_ef){
    t_max <- fu_max
    t <- 0
    s <- 0
    Lambda <- function(tupper) integrate(f = lambda,x=x,beta_x_rec=beta_x_rec, rnd_ef, lower = 0, upper = tupper)$value
    Lambda_inv <- function(s){
      v <- seq(0,t_max+3, length.out = 250)
      min(v[Vectorize(Lambda)(v)>=s])
    }

    X <- numeric(0)
    while(t <= t_max){
      u <- runif(1)
      s <- s -log(u)
      t <- Lambda_inv(s)
      X <- c( X, t)
    }

    return(X)
  }



  ## GERA TEMPO DE OCORRÊNCIA DOS EVENTOS ====
  #' @title gen_data_plp1
  #' @aliases gen_data_plp1
  #' @export
  #' @description           Gera os tempos de ocorrência e recorrência de todos os eventos.
  #' @param ID              Identificador dos indivíduos.
  #' @param N               Número de indivíduos.
  #' @param fu_max          tempo máximo de acompanhamento.
  #' @param x_cov           Matrix de covariáveis.
  #' @param rnd_ef_tot      Vetor de efeitos aleatórios.
  #' @param recurr          Indicador de recorrência.
  #'

  gen_data_pol_int<-function( ID,N, fu_max, x_cov,
                              beta_x_rec,
                              rnd_ef_tot,
                              recurr
  ){

    # Declara funćão lambda ----
    b <-10
    lambda_cov <- function(t,x,beta_x_rec,rnd_ef) (1 + b*(1+sin(0.25*pi*t)))*exp(t(x)%*%beta_x_rec + rnd_ef)

    # Gera tempos de ocorrências dos eventos ----

    tab_0<-c()
    for (i in 1:N){
      print(i)
      # x1<-X[i,2:3]
      x<-x_cov[i,]
      rnd_ef<-rnd_ef_tot[i]
      # print(x)
      if(recurr[i]==0){
        res_1 <-as.data.frame(get_nhpp_realization(lambda_cov,fu_max, x, beta_x_rec,rnd_ef))
      }else{
        res_1 <-as.data.frame(get_nhpp_realization_0(lambda_cov,fu_max, x, beta_x_rec,rnd_ef))
      }
      # print(head(res_1))
      # dim(res_1)
      colnames(res_1)<-c("end")
      res_1$ID<-i
      tab_0<-rbind(tab_0,res_1)
    }

    # cov<-cov.fu$x1

    # Organiza dados para entrada no spnhppzi ----
    tab<-tab_0%>%
      group_by(ID) %>%
      mutate(
        # seq=1:n(),
        # end=case_when(end>10~10,TRUE~end),
        ngroup=n(),
        rep=row_number(),
        expand=case_when((ngroup==rep & !(rep==1&(end==0 | end==fu_max)))~2,TRUE~1),
        expand1=expand)%>%
      expandRows("expand") %>%
      mutate(ngroup1=n(),
             IndRec=case_when(ngroup1>2~1, TRUE~0),
             rep1=row_number(),
             begin=case_when(rep1==1~0,TRUE~lag(end)),
             end=case_when(ngroup1==rep1~fu_max, TRUE~end),
             status=case_when(end==fu_max~0,TRUE~1))
    # %>%
    #   ungroup() %>%
    #   # left_join(cov, by="ID") %>%
    #   left_join(input_gen_data,by="ID")

    return(tab)
  }


