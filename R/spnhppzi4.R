######################################################################################################
#' @title spnhppzi
#' @aliases spnhppzi
#' @export
#' @description NHPPZI: Ajusta um modelo de processo de poisson não homegêneo
#              considerando a forma da Lei de Potência para função de intensidade e dados com
#              distribuição Inflada de Zeros. Permite abordagem frequentista e bayesiana,
#              além do uso de fragilidade compartilhada, onde rho(t) = rho0(t)*exp(X*beta + w)
#              ou rho(t) = rho0(t)*Z*exp(X*beta).
#' @param data: data frame contendo conjunto de dados
#' @param id: identificador de indivíduo
#' @param status: evento ou censura
#' @param stop: tempo até ocorrência do evento
#' @param IndRec: Indicador se o indivíduo apresenta uma ou mais recorrência
#' @param covar: Covariável
#' @param initial: Valores de incialização dos parâmetros
#' @param approach: abordagem a ser considerada (0 - frequentista, 1 - bayesiana)
#' @param frag:  uso de fragilidade (0 - sem fragilidade, 1 - com fragilidade)
#@tpfrag: tipo de fragilidade (0-rho(t) = rho0(t)*exp(X*beta + w), 1 - rho(t) = rho0(t)*Z*exp(X*beta))
#' @export
######################################################################################################
spnhppzi4<-function(formula,
                       data,
                       baseline = c("plp1", "plp2","plp3","bp"),
                       approach = c("mle", "bayes"),
                       n_iter=4000,
                       n_cores=1,
                       n_chains=4,
                       ZI = c("true","false"),
                       FR= c("true","false"),
                       initial,
                       tp_prior=0,frag=0,mu_omega=0,
                       sigma_omega=0,
                       #shp_sigma_omega=0, scl_sigma_omega=0,
                       shp_sigma2_z=0, scl_sigma2_z=0,
                       shp_alpha1=0,scl_alpha1=0,shp_alpha2=0,scl_alpha2=0,
                       mu_beta=0,sigma_beta=10,
                       mu_psi=0,sigma_psi=10,
                       rnd_logist=0,
                       mu_xi=0,sigma_xi=1,
                       spatial=0,
                       sp_model = c("car","sparse","icar"),
                       nb_mat=NULL,
                       shp_tau=0,
                       scl_tau=0,
                       W_n=0,
                       bp_degree=NULL,
                       h1_gamma=0,
                       h2_gamma=4
                       ){
  formula <- Formula::Formula(formula)
  approach <- tolower(approach)
  approach <- match.arg(approach)
  baseline <- tolower(baseline)
  baseline <- match.arg(baseline)
  ZI <- tolower(ZI)
  ZI <- match.arg(ZI)
  FR <- tolower(FR)
  FR <- match.arg(FR)
  sp_model<-tolower(sp_model)
  sp_model<-match.arg(sp_model)
  mf <- stats::model.frame(formula=formula, data=data)
  Terms <- stats::terms(mf)
  resp <- stats::model.response(mf)
  X <- stats::model.matrix(formula, data = mf, rhs = 1)[,-1, drop = FALSE]
  Z <- stats::model.matrix(formula, data = mf, rhs = 2)[,-1, drop = FALSE]
  time <- resp[,1]
  event <- resp[,2]
  id <- resp[,3]
  SP_ID<-resp[,4]
  IndRec<- resp[,5]


  data1<-as.data.frame(cbind(time,event,id,IndRec,SP_ID,X)) %>%
    group_by(id) %>%
    summarise(across(everything(), last,.names = "max_{.col}"))


  data2<-as.data.frame(cbind(time,event,id,IndRec,Z)) %>%
    group_by(id) %>%
    summarise(across(everything(), last,.names = "max_{.col}"))

  n_ind1<-as.data.frame(cbind(event,id)) %>%
    group_by(id) %>%
    summarise(n_ind1=sum(event)) %>%
    dplyr::select("n_ind1")
  n_ind1<-as.vector(unlist(n_ind1))
  # print(head(n_ind1,n=200))
  # print(max(n_ind1))

  position_ind<-as.data.frame(id) %>%
    group_by(id) %>%
    summarise(n_ind=n()) %>%
    mutate(end_ind=cumsum(n_ind), begin_ind=end_ind-n_ind+1) %>%
    dplyr::select(n_ind,begin_ind,end_ind)

  n_ind<-as.vector(unlist(position_ind$n_ind))
  begin_ind<-as.vector(unlist(position_ind$begin_ind))
  end_ind<-as.vector(unlist(position_ind$end_ind))

  # print(head(n_ind))

  # n_ind1<-as.data.frame(id) %>%
  #   group_by(id) %>%
  #   summarise(n_ind1=sum(event)) %>%
  #   select(n_ind1)
  # n_ind1<-as.vector(unlist(n_ind1))
  # print(head(n_ind1))

  max_stop <- as.vector(unlist(data1[,2]))
  IndRec2<- as.vector(unlist(data1[,4]))
  gr_SP_ID<- as.vector(unlist(data1[,5]))
  zeta<-max(time)

  N <- length(time)
  SP_N<-nrow(nb_mat)
  n <- length(max_stop)
  # print(length(max_stop))
  p <- ncol(X)
  q <- ncol(Z)
  print(q)

  #print(p)
  #print(q)

  if(p==0){
    X <- array(0, dim = c(0, 0))
    Xy <- array(0, dim = c(0, 0))
  }else{
    Xy<-as.matrix((data1[,-(1:5)] ))
  }
  #print(head(X))
  #print(head(Xy))
  if(q==0){
    Z <- array(0, dim = c(0, 0))
    Z1 <- array(0, dim = c(0, 0))
  }else{
    Z1<-as.matrix((data2[,-(1:4)] ))
  }
  #print(head(Z))
  # print(head(Z1,n=20L))
  baseline <- switch(baseline,
                     "plp" = 1,
                     "plp1" = 2,
                     "plp2" = 3,
                     "bp" = 4
  )

  if(baseline == 1 | baseline == 2|baseline == 3){
    m <- 2
  }

  g<-NULL
  G<-NULL
  if(baseline == 4){
    # zeta<-max(time)
    # print(zeta)
    bases <- bp(time, max_stop, bp_degree, zeta, N, n)
    g <- bases$b
    G <- bases$B
    m <- bp_degree
    # print(head(g))
    # print(head(G))
  }

  approach <- switch(approach,
                     "mle" = 0,
                     "bayes" = 1
  )

  FR <- switch(FR,
               "true" = 1,
               "false" = 0
  )

  ZI <- switch(ZI,
               "true" = 1,
               "false" = 0
  )

  sp_model<-switch(sp_model,
                   "car" = 0,
                   "sparse" = 1,
                   "icar" = 2
 )

  data_model <- list(id=id,evento=event,time=time, X=X, Z=Z, N=N, Xy=Xy, Z1=Z1,
                     max_stop=max_stop, n=n, p=p, q=q, IndR=IndRec, IndRec2=IndRec2, approach=approach, FR=FR, ZI=ZI,
                     begin_ind=begin_ind,end_ind=end_ind,
                     n_ind=n_ind,n_ind1=n_ind1, m=m, mu_omega=mu_omega,shp_sigma2_z=shp_sigma2_z, scl_sigma2_z=scl_sigma2_z,
                     sigma_omega=sigma_omega,
                     shp_alpha1=shp_alpha1, scl_alpha1=scl_alpha1, shp_alpha2=shp_alpha2, scl_alpha2=scl_alpha2,mu_beta=mu_beta,
                     sigma_beta=sigma_beta,mu_psi=mu_psi,
                     sigma_psi=sigma_psi, baseline=baseline, tp_prior=tp_prior,
                     SP_ID,gr_SP_ID, SP_N, nb_mat=nb_mat, shp_tau=shp_tau, scl_tau=scl_tau, W_n=W_n,
                     rnd_logist=rnd_logist,
                     mu_xi=mu_xi, sigma_xii=sigma_xi,
                     G=G, g=g, zeta=zeta
                     )
  if(spatial==0){
  if(FR==0){
    if(baseline==4){
      mod <-stanmodels$BPNHPP_COV_4

    }
    else{
    if(ZI==0){
      mod <- stanmodels$NHPP_COV_4
    }
    else{
      if(q==0){
        mod<- stanmodels$NHPP_ZI_03_07_2022
      }
      else{
        mod<- stanmodels$NHPP_ZI_LOGISTCOV_23_07_2022
      }
    }
    }
  }
  else {
    if(ZI==0){
  # mod <- rstan::stan_model("~/R/x86_64-pc-linux-gnu-library/3.6/spnhppzi/stan/NHPP_COV_FRAT_10.stan")
      #
      #mod <- rstan::stan_model("~/R/x86_64-pc-linux-gnu-library/4.2/NHPPZISP/stan/NHPP_COV_FRAT_10.stan")
      #mod<- rstan::stan_model("/usr/local/lib/R/site-library/NHPPZISP/stan/NHPP_COV_FRAT_10.stan")
      mod <- stanmodels$NHPP_COV_FRAT_10
    }
    else{
      if(q==0){
       #mod <- rstan::stan_model("~/R/x86_64-pc-linux-gnu-library/3.6/spnhppzi/stan/NHPP_ZI_FRAT_04_07_2022.stan")
      #
        #mod<- rstan::stan_model("~/R/x86_64-pc-linux-gnu-library/4.2/NHPPZISP/stan/NHPP_ZI_FRAT_09_03_2022.stan")
       # mod<- rstan::stan_model("/usr/local/lib/R/site-library/NHPPZISP/stan/NHPP_ZI_FRAT_09_03_2022.stan")
       #mod<- stanmodels$NHPP_ZI_FRAT_09_03_2022
        #mod<- rstan::stan_model("/home/alisson/R/x86_64-pc-linux-gnu-library/4.2/spnhppzi/stan/NHPP_ZI_FRAT_04_07_2022.stan")
       # mod<- rstan::stan_model("/home/alisson/spnhppzi/inst/stan/NHPP_ZI_FRAT_04_07_2022_1.stan")
        mod<- stanmodels$NHPP_ZI_FRAT_04_07_2022
      }
      else{
      #mod <- rstan::stan_model("~/R/x86_64-pc-linux-gnu-library/3.6/spnhppzi/stan/NHPP_ZI_LOGISTCOV_FRAT_16_03_2022.stan")
               #mod<- rstan::stan_model("~/R/x86_64-pc-linux-gnu-library/4.2/NHPPZISP/stan/NHPP_ZI_LOGISTCOV_FRAT_16_03_2022.stan")
        #mod<- rstan::stan_model("/usr/local/lib/R/site-library/NHPPZISP/stan/NHPP_ZI_LOGISTCOV_FRAT_16_03_2022.stan")
        mod<- stanmodels$NHPP_ZI_LOGISTCOV_FRAT_31_07_2022
      }
    }
  }
  }
  else{
    if(ZI==0){
    if(sp_model==0){
     print("CAR model")
     mod <- stanmodels$SPNHPP_COV_FRAT
     }
    if(sp_model==1){
      print("CAR Sparse model")
      mod <- stanmodels$SPNHPP_SPARSE_COV_FRAT
    }
    if(sp_model==2){
     print("ICAR model")
     mod <- stanmodels$SPNHPP_ICAR_COV_FRAT_15_09
    }
   # mod <- stanmodels$SPNHPP_COV_FRAT_10
    }
    else{
      if(q==0){
      print("ICAR model ZI")
      # mod <- stanmodels$SPNHPP_ZI_FRAT_15_09_2022
        mod <- stanmodels$SPNHPP_ZI_FRAT_05_11_2022
      }
      else{
      print("ICAR model ZI logistcov")
      mod <- stanmodels$SPNHPP_ZI_LOGISTCOV_FRAT_EFECT_IN_LOG_12_10_2022

      }
    }
  }
  if(FR==0){
    #Optimizing
    if(approach==0){
      result<- optimizing(stanmodels$mod, data = data_model, hessian =TRUE, init=initial, algorithm = "LBFGS")
      return(result)
    }

    if(approach==1){
      #result_b<- sampling(mod, data = data_model, cores = 4, iter=4000)
      result_b<- sampling(mod, data = data_model, cores = n_cores, iter=n_iter, chains=n_chains)
      # result_b<- sampling(mod, data = data_model, cores = 4, iter=4000,  control = list(max_treedepth = 50,adapt_delta = 0.999))
      return(result_b)
    }
  }
  else{
    #result_c<- sampling(mod, data = data_model, cores = 4, iter=4000)
   result_c<- sampling(mod, data = data_model, cores = n_cores, iter=n_iter,chains=n_chains)
    #result_c<- sampling(mod, data = data_model, cores = 4, iter=4000, control = list(max_treedepth = 50,adapt_delta = 0.999))
    return(result_c)
  }
}
