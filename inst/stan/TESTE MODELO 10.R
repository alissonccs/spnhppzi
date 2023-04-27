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
library(spam)
library(tibble)
library(loo)


folder="/home/alisson/UFMG/TESE/SIM/SIMULACOES MAR 2023/PC_ALISSON/Modelo_09_SP_BP/"
model_name<-"MOD10_SP"



model<-"Model 10 ICAR BP"
geo_level<-"Municipio + AP BH"
model_specif<-"Adj equal Gen"
N<-1000
alpha1_r<-2
alpha2_r<-1.3
# beta1_r<-0
# beta2_r<-0
beta1_r<-0.6
beta2_r<-0.8
sp_tau_r<-1
#sp_alpha_r<-0.99
psi1_r<-1.6
# psi2_r<-1.2
pi_r<-0.75
fu.min<-7
fu.max<-7
degree_bp<-min(ceiling(N^0.4),5)
# degree_bp<-min(20)
n_points_graf<-50
lista<-c(1:300)



set.seed(5832)
cov.fu<-spnhppzi::gencovfu2(N=N,
                            fu.min=fu.min,
                            fu.max=fu.max,
                            cens.prob = 0,
                            dist.x = c("binomial","normal","normal"),
                            par.x=list(0.7, c(0, 1),c(0, 2)),
                            beta.x=c(beta1_r,beta2_r,psi1_r)
)
set.seed(NULL)



base_sp<-spnhppzi::spsimrec8.1(N=cov.fu$N,
                     cov_rec=c("ID","X1","X2"),
                     cov_log=c("X3"),
                     beta_x_rec = c(beta1_r,beta2_r),
                     beta_x_log = c(psi1_r),
                     logist = 1,
                     x1 = cov.fu$x1,
                     fu = cov.fu$fu,
                     fu_max = cov.fu$fu.max,
                     fu_min = cov.fu$fu.min,
                     spatial=1,
                     list_area=list_area_RMBH,
                     sp_model="ICAR",
                     SP_N=100,
                     nb_mat=RMBH_mat,
                     sp_tau=sp_tau_r,
                     #sp_alpha=sp_alpha_r,
                     #beta.x = cov.fu$beta.x,
                     # dist_z = "lognormal",
                     random_ef=1,
                     # tp_rnd_ef=0,
                     pi = pi_r,
                     par_z=0,
                     dist_int_func = "weibull",
                     par_int_func = c(alpha1_r,alpha2_r),
                     baseline="plp2"
)



formula2=Formula(spnhppzi::Recur1(end,status,ID,SP_ID,IndRec)~X1+X2|X3)
RESULT_BAYES_SCOV1<- spnhppzi_model_04(formula2,
                                         base_sp,
                                         baseline = "plp2",
                                         FR="TRUE",
                                         ZI="TRUE",
                                         approach = "BAYES",
                                         sp_model = "ICAR",
                                         initial=1,
                                         tp_prior=1,
                                         shp_alpha1=1,scl_alpha1=1,
                                         shp_alpha2=1,scl_alpha2=1,
                                         mu_beta=0,sigma_beta=1,
                                         mu_psi=0,sigma_psi=1,
                                         mu_omega=0,
                                         # shp_sigma2_z =0.1, scl_sigma2_z =0.1,
                                         spatial=1,
                                         nb_mat=RMBH_mat,
                                         shp_tau=1,
                                         scl_tau=1,
                                         n_iter = 2000,
                                         n_cores=2,
                                         n_chains=2,
                                         W_n=284,omega_data = 0,
                                         # bp_degree=degree_bp,
                                         # h1_gamma=0,
                                         # h2_gamma=1,
                                         lower_tau=0,
                                         tp_prior_tau="gamma",
                                         tp_icar = 1,
                                         std_dev = 1
                                         # omega=omega_data$omega,
                                         # tau=sp_tau_r,
)



set.seed(NULL)
# pars<- as.data.frame(RESULT_BAYES_SCOV1$result_stan, pars = c("alpha","beta","pii","tau"))
pars_desc<-summary(RESULT_BAYES_SCOV1$result_stan, pars = c("alpha","beta","tau","psi"))
pars_desc<-pars_desc$summary
loglik<-as.data.frame(RESULT_BAYES_SCOV1$result_stan, pars = c("log_lik"))

omega_desc<-summary(RESULT_BAYES_SCOV1, pars = c("omega"))
omega_desc<-omega_desc$summary

tau_desc<-summary(RESULT_BAYES_SCOV1, pars = c("tau"))
tau_desc<-tau_desc$summary

# gamma_vec<-summary(RESULT_BAYES_SCOV1$result_stan, pars = c("gamma"))
# gamma_vec<-gamma_vec$summary
# gamma_vec<-gamma_vec[,1]

## CRITÉRIOS  ----
# WAIC
loo_WAIC<-loo::waic((as.matrix(loglik)))

loo_WAIC_estimates<-loo_WAIC$estimates %>%
  as.data.frame() %>%
  mutate(step=step)

# PSIS
loo_PSIS<-loo(RESULT_BAYES_SCOV1$result_stan,
              pars = "log_lik",
              save_psis = FALSE,
              cores = getOption("mc.cores", 1),
              moment_match = FALSE,
              k_threshold = 0.7)




###################################

## TRATA DADOS GEOGRÁFICOS ----
## Lista de códigos dos municípios da região metropolitana ====
list_RMBH<-c(3165537,3141108,3126000,3140704,3137601,3136652,3132206,3130101,
             3129806,3157807,3155306,3154804,3154606,3153905,3124104,3171204,
             3168309,3162955,3156700,3133709,3134608,3136603,3140159,3144805,
             3106705,3110004,3149309,3105004,3106200,3109006,3112505,3117876,
             3118601, 3162922
             # ,
             # #COLAR METROPOLITANO
             # 3157203, 3126406, 3142304, 3161908, 3108107, 3127206, 3106408, 3147105,
             # 3131901, 3167202, 3153608, 3107703, 3133808, 3163102, 3105400, 3131000
)


## SHAPE REGIÃO METROPOLITANA DE BH ====
MUN_MG<-st_read('/home/alisson/UFMG/TESE/SIM/MALHAS/MG_Municipios_2021/MG_Municipios_2021.shp')
RMBH<-MUN_MG[MUN_MG$CD_MUN %in% list_RMBH,]
dim(RMBH)
# mapview(RMBH)


### Exclui BH ====
RMBH1<-RMBH %>%
  filter(NM_MUN!="Belo Horizonte") %>%
  dplyr::select(!c(SIGLA,AREA_KM2))
nrow(RMBH1)

#### Shape áreas de ponderação de BH ====
APBH<-st_read('/home/alisson/UFMG/TESE/SIM/MALHAS/31_MG_Minas_Gerais/BELO HORIZONTE_area de ponderacao.shp')

APBH<-APBH %>%
  rename(CD_MUN=CD_APONDE) %>%
  mutate(NM_MUN="Belo Horizonte") %>%
  dplyr::select(!c(Qt_Setores))

# %>%
#   as.data.frame()

### Adiciona áreas de ponderação de BH ao shape da Região Metropolitana ====
RMBH2<-rbind(RMBH1,APBH)

RMBH2_sorted <- RMBH2[ order(as.numeric(RMBH2$CD_MUN)), ]
RMBH2_sorted$SP_ID<-seq(1:nrow(RMBH2_sorted))
RMBH_mat <- nb2mat(poly2nb(RMBH2_sorted), style = "B")
row.names(RMBH_mat)<-RMBH2_sorted$SP_ID
sum(RMBH_mat)/2
list_area_RMBH<-as.numeric(row.names(RMBH_mat))
