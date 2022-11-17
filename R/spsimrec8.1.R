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
spsimrec8.1 <-  function(N,
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
  source("/home/alisson/spnhppzi/R/Utils_spsimrec2.R")
  ID <- c(1:N)
  dist_z <- tolower(dist_z)
  dist_z <- match.arg(dist_z)
  nr.cov_rec <- length(beta_x_rec)
  print(nr.cov_rec)
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
                  pi,
                  logist,
                  cov_log,
                  beta_x_log,
                  xi,
                  rnd_ef1$rnd_ef)
   pi<-recurr$pi
   #print(pi)
   recurr<-recurr$recurr
   #print(recurr)
   recurr1<-as.data.frame((cbind(ID,recurr,pi)))
   colnames(recurr1)<-c("ID","recurr","pi")
   recurr_out<-list(recurr=recurr, recurr1=recurr1)


print(head(cov_rec))
print(head(cov_log))
   input_gen_data<-x1 %>%
     # left_join(cov_log,by="ID") %>%
     left_join(rnd_ef_out$rnd_ef1,by="ID") %>%
     left_join(recurr_out$recurr1,by="ID") %>%
     left_join(fu,by="ID")
   # %>%
   #   left_join(x1,by="ID")

   # print(head(fu))
  print(head(input_gen_data))
   # print(fu.max)
   print(head(input_gen_data$rnd_ef))
   print(head(rnd_ef1))
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
                  nr.cov_rec=nr.cov_rec)
  #tab <-gen_data(ID, N, dist_int_func, par_int_func, fu, x,rnd_ef)
   print(head(tab))
  base<-tab %>%
    left_join(input_gen_data,by="ID")

  return(base)

}
