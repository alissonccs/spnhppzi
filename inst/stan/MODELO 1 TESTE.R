

#######################################
#######################################

# remove.packages("spnhppzi")
# devtools::install_github("alissonccs/spnhppzi@master", auth_token = "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")

# install.packages("devtools")
library(devtools)
library(spnhppzi)
library(rstan)
library(Formula)
library(dplyr)
library(psych)
library(splitstackshape)
library(sf)
library(mapview)
library(spdep)
library(MASS)
library(matlib)
library(snowfall)
# install.packages("spam")
library(spam)
library(tibble)
model<-"Model 01"
# geo_level<-"Municipio + AP BH"
# model_specif<-"Adj equal Gen"
N<-1000
alpha1_r<-2
alpha2_r<-1.3
beta1_r<-0
beta2_r<-0
# beta1_r<-0.6
# beta2_r<-0.8
# sp_tau_r<-1
#sp_alpha_r<-0.99
# psi1_r<-1.6
# psi2_r<-1.2
pi_r<-0
fu.min<-7
fu.max<-7
# degree_bp<-min(ceiling(N^0.4),20)
# n_points_graf<-50

## Lista de códigos dos municípios da região metropolitana ====

# list_RMBH<-c(3165537,3141108,3126000,3140704,3137601,3136652,3132206,3130101,
#              3129806,3157807,3155306,3154804,3154606,3153905,3124104,3171204,
#              3168309,3162955,3156700,3133709,3134608,3136603,3140159,3144805,
#              3106705,3110004,3149309,3105004,3106200,3109006,3112505,3117876,
#              3118601, 3162922
#              # ,
#              # #COLAR METROPOLITANO
#              # 3157203, 3126406, 3142304, 3161908, 3108107, 3127206, 3106408, 3147105,
#              # 3131901, 3167202, 3153608, 3107703, 3133808, 3163102, 3105400, 3131000
# )
# 
# # length(list_RMBH)
# #RMBH<-munic_MG[munic_MG$CD_MUN %in% list_RMBH,]
# # plot(RMBH$geometry)
# # RMBH$geometry
# 
# ## SHAPE REGIÃO METROPOLITANA DE BH ====
# MUN_MG<-st_read('/home/alisson/UFMG/TESE/SIM/MALHAS/MG_Municipios_2021/MG_Municipios_2021.shp')
# RMBH<-MUN_MG[MUN_MG$CD_MUN %in% list_RMBH,]
# dim(RMBH)
# # mapview(RMBH)
# 
# 
# ### Exclui BH ====
# RMBH1<-RMBH %>%
#   filter(NM_MUN!="Belo Horizonte") %>%
#   dplyr::select(!c(SIGLA,AREA_KM2))
# nrow(RMBH1)
# 
# #### Shape áreas de ponderação de BH ====
# APBH<-st_read('/home/alisson/UFMG/TESE/SIM/MALHAS/31_MG_Minas_Gerais/BELO HORIZONTE_area de ponderacao.shp')
# 
# APBH<-APBH %>%
#   rename(CD_MUN=CD_APONDE) %>%
#   mutate(NM_MUN="Belo Horizonte") %>%
#   dplyr::select(!c(Qt_Setores))
# 
# nrow(APBH)
# # %>%
# #   as.data.frame()
# 
# ### Adiciona áreas de ponderação de BH ao shape da Região Metropolitana ====
# RMBH2<-rbind(RMBH1,APBH)
# row.names(RMBH2)
# # View(RMBH3)
# RMBH2_sorted <- RMBH2[ order(as.numeric(RMBH2$CD_MUN)), ]
# RMBH2_sorted$SP_ID<-seq(1:nrow(RMBH2_sorted))
# # View(RMBH3_sorted)
# # nrow(RMBH3)
# # list(RMBH3)
# RMBH_mat <- nb2mat(poly2nb(RMBH2_sorted), style = "B")
# row.names(RMBH_mat)<-RMBH2_sorted$SP_ID
# #View(RMBH_mat1)
# sum(RMBH_mat)/2
# list_area_RMBH<-as.numeric(row.names(RMBH_mat))
# 


# GERA COVARIÁVEIS ====
set.seed(5832)
cov.fu<-gencovfu2(N=N,
                  fu.min=fu.min,
                  fu.max=fu.max,
                  cens.prob = 0,
                  dist.x = c("binomial","normal"),
                  par.x=list(0, c(0, 0)),
                  beta.x=c(beta1_r, beta2_r)
)
set.seed(NULL)

# GERA DADOS ====

run_montecarlo <- function(step){
  set.seed(258258*step)
  
  base<-spsimrec8.1(N=cov.fu$N,
                    cov_rec=c("ID","X1","X2"),
                    # cov_log=c("X3","X4"),
                    beta_x_rec = c(beta1_r,beta2_r),
                    # beta_x_log = c(psi1_r,psi2_r),
                    logist = 0,
                    x1 = cov.fu$x1,
                    fu = cov.fu$fu,
                    fu_max = cov.fu$fu.max,
                    fu_min = cov.fu$fu.min,
                    spatial=0,
                    # list_area=list_area_RMBH,
                    # sp_model="ICAR",
                    # SP_N=100,
                    # nb_mat=RMBH_mat,
                    # sp_tau=sp_tau_r,
                    #sp_alpha=sp_alpha_r,
                    #beta.x = cov.fu$beta.x,
                    # dist_z = "lognormal",
                    random_ef=0,
                    tp_rnd_ef=0,
                    pi = pi_r,
                    par_z=0,
                    dist_int_func = "weibull",
                    par_int_func = c(alpha1_r,alpha2_r),
                    baseline="plp2"
                    
                    
  )
  # base_sp1<-base_sp %>%
  #   filter(recurr==0) %>%
  #   group_by(ID) %>%
  #   summarise(ngroup1=max(ngroup1),rnd_ef=max(rnd_ef),X1=max(X1),X2=max(X2))
  # 
  # 
  # media_sp<-base_sp1  %>%
  #   summarise(
  #     media = mean(ngroup1),
  #     mediana = median(ngroup1)
  #   )
  
  
  formula2=Formula(spnhppzi::Recur1(time=end,event=status,id=ID,SP_ID=NULL,IndRec=IndRec)~1|-1)
  RESULT_BAYES_SCOV1<- spnhppzi::spnhppzi4(formula2,
                                           base,
                                           baseline = "plp2",
                                           FR="FALSE",
                                           ZI="FALSE",
                                           approach = "BAYES",
                                           sp_model = "ICAR",
                                           initial=1,
                                           tp_prior=1,
                                           shp_alpha1=1,scl_alpha1=1,
                                           shp_alpha2=1,scl_alpha2=1,
                                           mu_beta=0,sigma_beta=0.5,
                                           mu_psi=0,sigma_psi=0.5,
                                           mu_omega=0,
                                           # shp_sigma2_z =0.1, scl_sigma2_z =0.1,
                                           spatial=0,
                                           # nb_mat=RMBH_mat,
                                           # shp_tau=1,
                                           # scl_tau=1,
                                           n_iter = 2000,
                                           n_cores=2,
                                           n_chains=2,
                                           # W_n=284,
                                           omega_data = 0,
                                           # bp_degree=degree_bp,
                                           # h1_gamma=0,
                                           # h2_gamma=1,
                                           # lower_tau=0,
                                           # tp_prior_tau="gamma",
                                           # tp_icar = 1,
                                           # std_dev = 1
                                           # omega=omega_data$omega,
                                           # tau=sp_tau_r,
                                           
                                           
  )
  
  set.seed(NULL)
  # pars<- as.data.frame(RESULT_BAYES_SCOV1, pars = c("alpha"))
  pars_desc<-summary(RESULT_BAYES_SCOV1,pars=c("alpha"))
  pars_desc<-pars_desc$summary
  # omega<-as.data.frame(RESULT_BAYES_SCOV1$result_stan, pars = c("omega"))
  # omega_des<-describe(omega,quant=c(.025,.5,.975),fast=TRUE)
  # tau<-as.data.frame(RESULT_BAYES_SCOV1$result_stan, pars = c("tau"))
  # tau_desc<-describe(tau,quant=c(.025,.5,.975),fast=TRUE)
  # 
  # 
  # gamma<-as.data.frame(RESULT_BAYES_SCOV1$result_stan, pars = c("gamma"))
  # gamma<-describe(gamma,quant=c(.025,.5,.975),fast=TRUE)
  # gamma_vec<-gamma[,3]
  # 
  # # 
  # max_time_id<-base_sp %>% 
  #   group_by(ID) %>% 
  #   summarise(max_time_id=max(end)) %>% 
  #   dplyr::select(max_time_id) %>% 
  #   as.vector()
  # max_time_id1<-as.vector(unlist(max_time_id))
  # # 
  # # 
  # # bp1 <- function(time, max_time_id, degree, zeta,N,n) {
  # #   # n <- length(time)
  # #   y <- time/zeta
  # #   y1 <-max_time_id/zeta
  # #   b <- matrix(nrow=N, ncol=degree)
  # #   B <- matrix(nrow=N, ncol=degree)
  # #   for(k in 1:degree)
  # #   {
  # #     # b[,k] <- stats::dbeta(y, k, degree - k + 1)/zeta
  # #     b[,k] <- stats::dbeta(y, k, degree - k + 1)
  # #     B[,k] <- stats::pbeta(y, k, degree - k + 1)
  # #   }
  # #   return(list(b=b, B=B))
  # # }
  # 
  # bases <- bp1(time=base_sp$end, max_time_id=max_time_id1, degree=degree_bp, zeta=7, N=length(base_sp$end), n=length(max_time_id1))
  # prod<-bases$B %*% (gamma_vec) 
  # 
  # data_bp<-cbind(prod,base_sp$end) %>% 
  #   as.data.frame() %>% 
  #   mutate(type="PB"
  #          ,
  #          step=0
  #          )
  # colnames(data_bp)<-c("FMA_P","time","type","step")
  # 
  # FMA<-(data_bp$time)^alpha2_r*alpha1_r
  # 
  # data_graf_PL<-cbind(FMA,data_bp$time) %>% 
  #   as.data.frame() %>% 
  #   mutate(type="PL",step=0)
  # 
  # colnames(data_graf_PL)<-c("FMA_P","time","type","step")
  # 
  # 
  # data_graf_total<-rbind(data_graf_PL,data_bp) 
  # 
  # # ggplot(economics, aes(x=date)) + 
  # #   geom_line(aes(y = psavert), color = "darkred") + 
  # #   geom_line(aes(y = uempmed), color="steelblue", linetype="twodash") 
  # 
  # ggplot2::ggplot()+
  #   geom_line(data=data_graf_total[data_graf_total$type!="PL",],aes(x=time,y=FMA_P,group=as.character(step)),colour="grey60")+
  #   geom_line(data=data_graf_total[data_graf_total$type=="PL",],aes(x=time,y=FMA_P),colour = "red")
  
  ############################################
  ############################################
  
  
  # table(data_graf_total$step)
  
  
  pars_desc <-pars_desc %>%
    as.data.frame %>%
    tibble::rownames_to_column("pars") %>%
    mutate(pars=case_when(pars=="alpha[1]"~"alpha1",
                          pars=="alpha[2]"~"alpha2"
                          # pars=="beta[1]"~"beta1",
                          # pars=="beta[2]"~"beta2",
                          # pars=="tau"~"sp_tau",
                          # pars=="pii[1]"~"pi",
                          # pars=="gamma[1]"~"gamma1",
                          # pars=="gamma[2]"~"gamma2",
                          # pars=="gamma[3]"~"gamma3",
                          # pars=="gamma[4]"~"gamma4",
                          # pars=="gamma[5]"~"gamma5"
                          
                          # ,
                          # pars=="sp_tau_r"~sp_tau_r,
    ),
    real_pars=case_when(pars=="alpha1"~alpha1_r,
                        pars=="alpha2"~alpha2_r
                        # pars=="beta1"~beta1_r,
                        # pars=="beta2"~beta2_r,
                        # pars=="sp_tau"~sp_tau_r,
                        # pars=="pi"~pi_r
                        # ,
                        # pars=="sp_tau_r"~sp_tau_r,
    ),
    step=step,
    model=model,
    # geo_level=geo_level,
    # model_specif=model_specif,
    sample=paste("Sample="," ",N,sep=""),
    period=7,
    # exp_ev<-(period)^alpha2_r*alpha1_r,
    exp_ev=paste("Expected events= ",round((period)^alpha2_r*alpha1_r),sep=""),
    bias= mean-real_pars,
    eqm= bias^2,
    #rbias=case_when(vars!=6~100*bias/abs(real_pars),TRUE~100*bias),
    rbias=100*bias/abs(real_pars),
    mod=paste("Model",model,", ZI=",pi_r,"%",", sample",sample, ", period=", period," years,", " E[N(t)]=",exp_ev))
  # pars_desc$step<-step
  # pars_desc$alpha1_r<-alpha1_r
  # pars_desc$alpha2_r<-alpha2_r
  # pars_desc$beta1_r<-beta1_r
  # pars_desc$beta2_r<-beta2_r
  # pars_desc$sp_tau_r<-sp_tau_r
  write.table(pars_desc, file="/home/alisson/UFMG/TESE/SIM/RESULT_25_03_2023/mod01_2.0_1.3_1m_12_03_2023.txt",
              row.names =TRUE,
              col.names =!file.exists("/home/alisson/UFMG/TESE/SIM/RESULT_25_03_2023/mod01_2.0_1.3_1m_12_03_2023.txt"),
              append = TRUE)
  
  # time_graf=seq(0,fu.max,fu.max/n_points_graf)
  # bases <- bp1(time=time_graf, max_time_id=rep(fu.max,length(time_graf)), degree=degree_bp, zeta=fu.max, N=length(time_graf), n=N)
  # prod<-bases$B %*% (gamma_vec)
  # 
  # data_bp<-cbind(prod,time_graf) %>%
  #   as.data.frame() %>%
  #   mutate(type="PB",
  #          step=step)
  # colnames(data_bp)<-c("FMA_P","time","type","step")
  # 
  # 
  # write.table(data_bp, file="/home/leste/Alisson/data_mod01_2.0_1.3_300_12_03_2023.txt",
  #             row.names =FALSE,
  #             col.names =!file.exists("/home/leste/Alisson/data_mod01_2.0_1.3_300_12_03_2023.txt"),
  #             append = TRUE)
  return(pars_desc)
}




cpus <- 4                       # UFMG
# cpus <- 1

sfInit(parallel=TRUE, cpus=cpus,slaveOutfile="/home/alisson/UFMG/TESE/SIM/RESULT_25_03_2023/log_mod01_2.0_1.3_1m_12_03_2023.txt")
sfExportAll()
sfLibrary(spnhppzi)
# sfLibrary(NHPZZISP::simrecev)
# sfLibrary(NHPZZISP::func_sim)
sfLibrary(dplyr)
sfLibrary(Formula)
sfLibrary(rstan)
sfLibrary(psych)
sfLibrary(splitstackshape)
sfLibrary(MASS)
sfLibrary(spam)
sfLibrary(tibble)
sfLapply(1:300, fun=run_montecarlo)
sfStop()
