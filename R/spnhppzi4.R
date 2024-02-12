######################################################################################################
#' @title spnhppzi
#' @aliases spnhppzi
#' @export
#' @description NHPPZI:
#' @param data: data frame contendo conjunto de dados
#' @param id: identificador de indivíduo
#' @param status: evento ou censura
#' @param stop: tempo até ocorrência do evento
#' @param IndRec: Indicador se o indivíduo apresenta uma ou mais recorrência
#' @param covar: Covariável
#' @param initial: Valores de incialização dos parâmetros
#' @param approach: abordagem a ser considerada (0 - frequentista, 1 - bayesiana)
#' @param frag:  uso de fragilidade (0 - sem fragilidade, 1 - com fragilidade)
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
                       tp_prior=0,
                       frag=0,
                       tp_rnd_ef=0,
                       mu_omega=0,
                       sigma_omega=0,
                       shp_sigma2_z=0, scl_sigma2_z=0,
                       shp_alpha1=0,scl_alpha1=0,shp_alpha2=0,scl_alpha2=0,
                       mu_beta=0,sigma_beta=10,
                       mu_psi=0,sigma_psi=10,
                       rnd_logist=0,
                       mu_xi=0,sigma_xi=1,
                       spatial=0,
                       sp_model = c("car","sparse","icar"),
                       nb_mat=NULL,
                       data_tau=0,
                       tau=1, #usar somente para testar modelo onde tau é data
                       tp_prior_tau=c("gamma","inv_gamma"),
                       shp_tau=0,
                       scl_tau=0,
                       lower_tau=0,
                       tp_icar=0,
                       W_n=0,
                       bp_degree=NULL,
                       h1_gamma=0,
                       h2_gamma=4,
                       omega_data=0,
                       omega=NULL,
                       std_dev=1,
                       tp_DIC=0
                       ){

  formula <- Formula::Formula(formula)
  approach <- tolower(approach)
  approach <- match.arg(approach)
  baseline <- tolower(baseline)
  baseline <- match.arg(baseline)
  tp_prior_tau<-tolower(tp_prior_tau)
  tp_prior_tau<-match.arg(tp_prior_tau)
  tp_prior_tau <- switch(tp_prior_tau,
                         "gamma" = 1,
                         "inv_gamma" = 2)
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

  position_ind<-as.data.frame(id) %>%
    group_by(id) %>%
    summarise(n_ind=n()) %>%
    mutate(end_ind=cumsum(n_ind), begin_ind=end_ind-n_ind+1) %>%
    dplyr::select(n_ind,begin_ind,end_ind)

  n_ind<-as.vector(unlist(position_ind$n_ind))
  begin_ind<-as.vector(unlist(position_ind$begin_ind))
  end_ind<-as.vector(unlist(position_ind$end_ind))

  max_stop <- as.vector(unlist(data1[,2]))
  IndRec2<- as.vector(unlist(data1[,4]))
  gr_SP_ID<- as.vector(unlist(data1[,5]))
  zeta<-max(time)

  N <- length(time)
  SP_N<-nrow(nb_mat)
  n <- length(max_stop)
  p <- ncol(X)
  q <- ncol(Z)
  # print(q)

  if(p==0){
    X <- array(0, dim = c(0, 0))
    Xy <- array(0, dim = c(0, 0))
  }else{
    Xy<-as.matrix((data1[,-(1:5)] ))
  }

  if(q==0){
    Z <- array(0, dim = c(0, 0))
    Z1 <- array(0, dim = c(0, 0))
  }else{
    Z1<-as.matrix((data2[,-(1:4)] ))
  }

  baseline <- switch(baseline,
                     "plp1" = 1,
                     "plp2" = 2,
                     "plp3" = 3,
                     "bp" = 4
  )
  if(baseline == 1 | baseline == 2|baseline == 3){
    m <- 2
  }

  g<-NULL
  G<-NULL
  if(baseline == 4){
    bases <- bp(time, max_stop, bp_degree, zeta, N, n)
    g <- bases$b
    G <- bases$B
    m <- bp_degree
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

# Lista de argumentos para o rStan
  data_model <- list(id=id,evento=event,time=time, X=X, Z=Z, N=N, Xy=Xy, Z1=Z1,
                     max_stop=max_stop, n=n, p=p, q=q, IndR=IndRec, IndRec2=IndRec2,
                     approach=approach, FR=FR, ZI=ZI, begin_ind=begin_ind,end_ind=end_ind,
                     n_ind=n_ind,n_ind1=n_ind1, m=m, mu_omega=mu_omega, tp_rnd_ef=tp_rnd_ef,
                     shp_sigma2_z=shp_sigma2_z, scl_sigma2_z=scl_sigma2_z,sigma_omega=sigma_omega,
                     shp_alpha1=shp_alpha1, scl_alpha1=scl_alpha1, shp_alpha2=shp_alpha2,
                     scl_alpha2=scl_alpha2,mu_beta=mu_beta, sigma_beta=sigma_beta,mu_psi=mu_psi,
                     sigma_psi=sigma_psi, baseline=baseline, tp_prior=tp_prior,
                     SP_ID,gr_SP_ID, SP_N, nb_mat=nb_mat, shp_tau=shp_tau, scl_tau=scl_tau, W_n=W_n,
                     rnd_logist=rnd_logist, mu_xi=mu_xi, sigma_xii=sigma_xi, G=G, g=g, zeta=zeta,
                     omega=omega, tau=tau, lower_tau=lower_tau , tp_prior_tau=tp_prior_tau,tp_icar=tp_icar,
                     std_dev=std_dev
                     )
  # Modelos não espaciais ----
  if(spatial==0){
  ## Sem efeito aleatório para o indivíduo ----
  if(FR==0){
    ### Semiparamétrico ----
    if(baseline==4){
      ##### Sem inflação de zeros ----
      if(ZI==0){
      mod <-stanmodels$S_NHPP
      }
      #### Com inflaćão de zeros ----
      else{
          mod<- stanmodels$SZI_NHPP
      }
    }
   ### Paramétrico ----
   else{
    #### Sem inflaćão de zeros ----
    if(ZI==0){
      mod <- stanmodels$NHPP
    }
    #### Com inflaćão de zeros ----
    else{
        mod<- stanmodels$ZI_NHPP
    }
    }
  }
  ## Com efeito aleatório para o indivíduo ----
  else {
    #### Sem inflaćão de zeros ----
    if(ZI==0){
      #### Semiparamétrico ----
      if(baseline==4){
       mod <-stanmodels$S_NHPP_RE
      }
      #### Paramétrico ----
      else{
      if(tp_DIC==0){#Usa Transformed Parameters
      mod <- stanmodels$NHPP_RE
      }
      else{#Usa Generate Quantities
      mod <- stanmodels$NHPP_RE_1
      }
     }
    }
    #### Com inflaćão de zeros ----
    else{
      # Sem covariável na reg logística ---
      if(q==0){
          # Semiparamétrico
          if(baseline==4){
            mod <-stanmodels$SZI_NHPP_RE
          }
          # Paramétrico
          else{
           mod<- stanmodels$ZI_NHPP_RE
          }
      }
      # Com covariável na reg logística ----
      else{
        # Semiparamétrico ----
        if(baseline==4){
          mod<- stanmodels$SZI_NHPP_RE_COV
        }
        # Paramétrico ----
        else{
         mod<- stanmodels$ZI_NHPP_RE_COV
        }
      }
    }
  }
  }
  # Modelos não espaciais ----
  else{
    # Sem inflaćão de zeros ----
    if(ZI==0){
    if(sp_model==2){
     print("ICAR model")
    # Semiparamétrico ----
    if(baseline==4){
     mod <- stanmodels$S_NHPP_SE
    }
      # Paramétrico ----
      else{
       mod <- stanmodels$NHPP_SE
      }
    }
    }
    # Com inflaćão de zeros ----
    else{
      # Sem covariável na reg logística ----
      if(q==0){
       print("ICAR model ZI")
          # Semiparamétrico ----
          if(baseline==4){
           mod <- stanmodels$SZI_NHPP_SE
          }
          # Paramétrico ----
          else{
           mod <- stanmodels$ZI_NHPP_SE
          }
      }
      # Com covariável na reg logística
      else{
        # Semiparamétrico ----
        if(baseline==4){
          mod <- stanmodels$SZI_NHPP_SE_COV
        }
        # Paramétrico ----
        else{
         print("ICAR model ZI logistcov")
         mod <- stanmodels$ZI_NHPP_SE_COV
        }
      }
    }
  }
  if(FR==0){
    #Optimizing
    if(approach==0){
      result<- optimizing(stanmodels$mod, data = data_model, hessian =TRUE, init=initial, algorithm = "LBFGS")
      if(baseline==4){
       result<-list("result_stan"=result,"G"=G, "g"=g, "m"=m)
      }
      return(result)
    }
    if(approach==1){
     result_b<- sampling(mod, data = data_model, cores = n_cores, iter=n_iter, chains=n_chains)
      if(baseline==4){
       result_b<-list("result_stan"=result_b,"G"=G, "g"=g, "m"=m)
      }
      return(result_b)
    }
  }
  else{
   result_c<- sampling(mod, data = data_model, cores = n_cores, iter=n_iter,chains=n_chains)
   if(baseline==4){
    result_c<-list("result_stan"=result_c,"G"=G, "g"=g, "m"=m)
   }
    return(result_c)
  }
}
