#' @title spsimrec
#' @aliases spsimrec
#' @export
#' @description
#' Simulates recurrent event data with options for spatial correlation, zero inflation,
#' random effects, and different intensity functions. The function allows customization
#' of the follow-up period, covariates, and model parameters.
#'
#' @param N Integer. Number of individuals in the dataset.
#' @param spatial Logical (0 or 1). Indicator for spatial modeling (default = 0, non-spatial).
#' @param sp_model Character. Type of spatial model used (`"ICAR"` or `"CAR"`).
#' @param SP_N Integer. Number of area units when using a spatial model.
#' @param nb_mat Matrix. Neighborhood matrix defining spatial relationships.
#' @param sp_tau Numeric. Precision parameter (`τ`) for the ICAR model.
#' @param sp_alpha Numeric. Spatial association parameter for the CAR model.
#'
#' ## Covariates
#' @param x Matrix. Matrix of covariates.
#' @param x1 Matrix. Matrix containing `ID` and covariates.
#' @param cov_rec Character vector. List of covariates associated with the intensity function.
#' @param cov_log Character vector. List of covariates associated with the logistic regression.
#' @param beta_x_rec Numeric vector. Coefficients of covariates in the intensity function.
#' @param beta_x_log Numeric vector. Coefficients of covariates in the logistic regression.
#'
#' ## Follow-up and Zero Inflation
#' @param fu_min Numeric. Minimum follow-up time.
#' @param fu_max Numeric. Maximum follow-up time.
#' @param pi_zi Numeric. Proportion of individuals with zero recurrences (zero inflation).
#' @param cens.prob Numeric. Probability of censoring before the maximum follow-up period.
#'
#' ## Random Effects
#' @param random_ef Logical (0 or 1). Indicator for using random effects (default = 0, no random effects).
#' @param tp_rnd_ef Logical (0 or 1). Indicator for using multiplicative random effects.
#' @param dist_z Character. Distribution of the non-spatial random effect (`"Gamma"` or `"lognormal"`).
#' @param par_z Numeric vector. Variance parameters for the random effect distribution.
#'
#' ## Intensity Function
#' @param dist_int_func Character. Type of intensity function (`"Weibull"` for power-law process).
#' @param par_int_func Numeric vector. Scale and shape parameters of the intensity function.
#' @param xi Numeric. Parameter linking the intensity function and logistic regression.
#'
#' ## Logistic Regression Model
#' @param logist Logical (0 or 1). Indicator for using covariates in the logistic regression model.
#' @param dist.x Character. Distribution of covariates (`"binomial"` or `"normal"`).
#' @param par.x List. Parameters of the covariate distributions.

# SIMRECEV - SIMULAÇÃO DE EVENTOS RECORRENTES ====
spsimrec<-  function(N,
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
                         pi_zi = 0,
                         # dfree = 0,
                         logist = 0,
                         mu_omega=0,
                         sigma_omega=0,
                         baseline = c("plp","polynom")
){
  # source("/home/alisson/spnhppzi/R/Utils_spsimrec.R")
  ID <- c(1:N)
  dist_z <- tolower(dist_z)
  dist_z <- match.arg(dist_z)
  nr.cov_rec <- length(beta_x_rec)
  print(nr.cov_rec)
  nr.cov_log <- length(beta_x_log)
  cov_rec<-x1 %>% dplyr::select(all_of(cov_rec))
  cov_log<-x1 %>% dplyr::select(all_of(cov_log))

  baseline <- tolower(baseline)
  baseline <- match.arg(baseline)

  baseline <- switch(baseline,
                     "plp" = 1,
                     "polynom"=2
  )

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
  ### Executa funções para gerar efeitos aleatórios ====
  if(random_ef==0){  #sem efeitos
    rnd_ef<-rep(1,N)
    tp_rnd_ef<-0
    rnd_ef1<-as.data.frame(cbind(ID,rnd_ef))
    colnames(rnd_ef1)<-c("ID","rnd_ef")
    rnd_ef_out<-list(rnd_ef=rnd_ef,rnd_ef1=rnd_ef1)
  }else if(spatial==0){ #efeitos não espaciais
    rnd_ef<- gen_rnd_ef(N, ID, dist_z, tp_rnd_ef, par_z,mu_omega,sigma_omega)
    rnd_ef1<-as.data.frame(cbind(ID,rnd_ef))
    colnames(rnd_ef1)<-c("ID","rnd_ef")
    rnd_ef_out<-list(rnd_ef=rnd_ef,rnd_ef1=rnd_ef1)
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
    rnd_ef_out<-list(rnd_ef=rnd_ef,rnd_ef1=rnd_ef1)
    # rnd_ef<-rnd_ef$rnd_ef
  }

  # INFLAÇÃO DE ZEROS - Define indivíduos recorrentes  ====
  ### Executa função ====
  recurr<-gen_zi(ID,
                 N,
                 pi_zi,
                 logist,
                 cov_log,
                 beta_x_log,
                 xi,
                 rnd_ef1$rnd_ef)
  pi_zi<-recurr$pi_zi
  #print(pi_zi)
  recurr<-recurr$recurr
  #print(recurr)
  recurr1<-as.data.frame((cbind(ID,recurr,pi_zi)))
  colnames(recurr1)<-c("ID","recurr","pi_zi")
  recurr_out<-list(recurr=recurr, recurr1=recurr1)


  # print(head(cov_rec))
  # print(head(cov_log))
  input_gen_data<-x1 %>%
    # left_join(cov_log,by="ID") %>%
    left_join(rnd_ef_out$rnd_ef1,by="ID") %>%
    left_join(recurr_out$recurr1,by="ID") %>%
    left_join(fu,by="ID")

  if (baseline==2){

  tab<- gen_data_pol_int(ID=input_gen_data$ID,
                         N=N,
                         fu_max=fu_max,
                         x_cov=as.matrix(input_gen_data[,2:(1+nr.cov_rec)]),
                         beta_x_rec=beta_x_rec,
                         # tp_rnd_ef= tp_rnd_ef,
                         rnd_ef_tot=input_gen_data$rnd_ef,
                         recurr=input_gen_data$recurr
  )
  } else{

  ## TEMPO DE OCORRÊNCIA DOS EVENTOS ====
  ### Executa função ====
  tab <-gen_data(ID=input_gen_data$ID,N=N,
                 dist_int_func=dist_int_func,
                 par_int_func=par_int_func,
                 fu=input_gen_data$fu,
                 fu_max =fu_max,
                 x=as.matrix(input_gen_data[,2:(1+nr.cov_rec)]),
                 # x1=x1,
                 beta_x_rec=beta_x_rec,
                 tp_rnd_ef= tp_rnd_ef,
                 rnd_ef=input_gen_data$rnd_ef,
                 # rnd_ef1=rnd_ef1,
                 recurr=input_gen_data$recurr,
                 # recurr1 = recurr1,
                 nr.cov_rec=nr.cov_rec,
                 baseline=baseline)
  }
  #tab <-gen_data(ID, N, dist_int_func, par_int_func, fu, x,rnd_ef)
  # print(head(tab))
  base<-tab %>%
    left_join(input_gen_data,by="ID")

  return(base)

}
