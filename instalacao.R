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
#remove.packages("spnhppzi")


library("spnhppzi")
library(spnhppzi)
?spnhppzi::spsimrec
?spnhppzi::Recur
library(rstan)
library(Formula)
library(dplyr)
library(psych)
library(splitstackshape)


base1<-spnhppzi::spsimrec(
  N=1000,
  fu.min=5,
  fu.max=5,
  #dist.x=c("binomial", "normal"),
  #par.x=list(0.5, c(0, 1)),
  #beta.x=c(1,1.3),
  random.ef = 0,
  dist.z="lognormal",
  ent.dist.z=0,
  mu.omega=0,
  sigma.omega=1,
  par.z = 1,
  dist.rec="weibull",
  par.rec=c(1,1.3),
  cens.prob=0,
  dfree=0,
  pi=0,
  logist=0
)

base1$IndRec<-1-base1$recurr



formula2=Formula(Recur(end,status,ID,IndRec)~1|-1)
RESULT_BAYES_SCOV1<- spnhppzi::spnhppzi(formula2,
                                         base1,
                                         baseline = "plp1",
                                         FR="false",
                                         ZI="false",
                                         approach = "BAYES",
                                         initial=10,
                                         tp_prior=1,
                                         shp_alpha1=0.1,scl_alpha1=0.1,
                                         shp_alpha2=0.1,scl_alpha2=0.1,
                                         mu_beta=0,sigma_beta=1,
                                         mu_psi=0,sigma_psi=1,
                                         mu_omega=0,
                                         shp_sigma_omega = 1, scl_sigma_omega = 1
)
pars<- as.data.frame(RESULT_BAYES_SCOV1, pars = c("alpha"))
colnames(pars)<-c("alpha1","alpha2")
pars<-pars %>%
  mutate(alpha1_1=1/alpha1)
pars_desc<-describe(pars,quant=c(.025,.5,.975),fast=TRUE)
