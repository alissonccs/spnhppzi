# library("rstantools")
# rstan_create_package(path = '/home/alisson/spnhppzi')
# setwd("/home/alisson/spnhppzi")
# list.files(all.files = TRUE)
# file.show("DESCRIPTION")
# file.show("Read-and-delete-me")
# file.remove('Read-and-delete-me')
# file.show(file.path("R", "spnhppzi-package.R"))
example(source) # defines the sourceDir() function
try(roxygen2::roxygenize(load_code = sourceDir), silent = TRUE)
roxygen2::roxygenize()

install.packages("../spnhppzi", repos = NULL, type = "source")
library("spnhppzi")
?spnhppzi::spsimrec
?spnhppzi::spnhppzi
# remove.packages("spnhppzi")


library("spnhppzi")
library(spnhppzi)
?spnhppzi::spsimrec
?spnhppzi::Recur
library(rstan)
library(Formula)
library(dplyr)
library(psych)
library(splitstackshape)

# EXAMPLE ----

N<-500
alpha1_r<-0.5
alpha2_r<-1.3
beta1_r<-0.6
beta2_r<-0.8
pi_r<-0
fu.min<-7
fu.max<-7
lista<-c(1:10)

set.seed(5832)
cov.fu<-gencovfu2(N=N,
                  fu.min=fu.min,
                  fu.max=fu.max,
                  cens.prob = 0,
                  dist.x = c("binomial","normal"),
                  par.x=list(0.7, c(0, 1)),
                  beta.x=c(beta1_r, beta2_r)
)
set.seed(NULL)

base<-spsimrec(N=cov.fu$N,
                cov_rec=c("ID","X1","X2"),
                beta_x_rec = c(beta1_r,beta2_r),
                logist = 0,
                x1 = cov.fu$x1,
                fu = cov.fu$fu,
                fu_max = cov.fu$fu.max,
                fu_min = cov.fu$fu.min,
                spatial=0,
                random_ef=0,
                tp_rnd_ef=0,
                pi = pi_r,
                par_z=0,
                dist_int_func = "weibull",
                par_int_func = c(alpha1_r,alpha2_r),
                baseline="plp2"
)

formula2=Formula(spnhppzi::Recur1(time=end,event=status,id=ID,SP_ID=NULL,IndRec=IndRec)~X1+X2|-1)
RESULT_BAYES_SCOV1<- spnhppzi::spnhppzi(formula2,
                                         base,
                                         baseline = "plp2",
                                         FR="FALSE",
                                         ZI="FALSE",
                                         approach = "BAYES",
                                         sp_model = "ICAR",
                                         initial=1,
                                         tp_prior=1,
                                         shp_alpha1=0.1,scl_alpha1=0.1,
                                         shp_alpha2=0.1,scl_alpha2=0.1,
                                         mu_beta=0,sigma_beta=4,
                                         mu_psi=0,sigma_psi=4,
                                         mu_omega=0,
                                         spatial=0,
                                         n_iter = 2000,
                                         n_cores=2,
                                         n_chains=2,
                                         omega_data = 0
)


summary(RESULT_BAYES_SCOV1,pars=c("alpha","beta"))


