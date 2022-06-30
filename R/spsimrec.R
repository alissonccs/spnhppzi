#'@title spsimrec
#' @aliases spsimrec
#' @export
#' @description Recur:
#' @param cens.prob  Gives the probability of being censored due to loss to follow-up before
#'   \code{fu.max}. For a random set of individuals defined by a B(N,\code{cens.prob})-distribution,
#'   the time to censoring is generated from a uniform
#'   distribution on \code{[0, fu.max]}. Default is \code{cens.prob=0}, i.e. no censoring
#'   due to loss to follow-up.
#' @param dist.x     Distribution of the covariate(s) \eqn{X}. If there is more than one covariate,
#'   \code{dist.x} must be a vector of distributions with one entry for each covariate. Possible
#'   values are \code{"binomial"} and \code{"normal"}, default is \code{dist.x="binomial"}.
#' @param par.x      Parameters of the covariate distribution(s). For \code{"binomial", par.x} is
#'   the probability for \eqn{x=1}. For \code{"normal"}, \code{par.x=c(}\eqn{\mu, \sigma}\code{)}
#'   where \eqn{\mu} is the mean and \eqn{\sigma} is the standard deviation of a normal distribution.
#'   If one of the covariates is defined to be normally distributed, \code{par.x} must be a list,
#'   e.g. \code{ dist.x <- c("binomial", "normal")} and \code{par.x  <- list(0.5, c(1,2))}.
#'   Default is \code{par.x=0}, i.e. \eqn{x=0} for all individuals.
#' @param beta.x  Regression coefficient(s) for the covariate(s) \eqn{x}. If there is more than one
#'   covariate, \code{beta.x} must be a vector of coefficients with one entry for each covariate.
#'   \code{simrec} generates as many covariates as there are entries in \code{beta.x}. Default is
#'   \code{beta.x=0}, corresponding to no effect of the covariate \eqn{x}.
#' @param dist.z     Distribution of the frailty variable \eqn{Z} with \eqn{E(Z)=1} and
#'   \eqn{Var(Z)=\theta}. Possible values are \code{"gamma"} for a Gamma distributed frailty
#'    and \code{"lognormal"} for a lognormal distributed frailty.
#'    Default is \code{dist.z="gamma"}.
#' @param par.z      Parameter \eqn{\theta} for the frailty distribution: this parameter gives
#'   the variance of the frailty variable \eqn{Z}.
#'   Default is \code{par.z=0}, which causes \eqn{Z=1}, i.e. no frailty effect.
#' @param dist.rec   Form of the baseline hazard function. Possible values are \code{"weibull"} or
#'   \code{"gompertz"} or \code{"lognormal"} or \code{"step"}.
#' @param par.rec  Parameters for the distribution of the event data.
#'   If \code{dist.rec="weibull"} the  hazard function is \deqn{\lambda_0(t)=\lambda*\nu* t^{\nu - 1},}
#'   where \eqn{\lambda>0} is the scale and \eqn{\nu>0} is the shape parameter. Then
#'   \code{par.rec=c(}\eqn{\lambda, \nu}\code{)}. A special case
#'   of this is the exponential distribution for \eqn{\nu=1}.\\
#'   If \code{dist.rec="gompertz"}, the hazard function is \deqn{\lambda_0(t)=\lambda*exp(\alpha t),}
#'   where \eqn{\lambda>0} is the scale and \eqn{\alpha\in(-\infty,+\infty)} is the shape parameter.
#'   Then \code{par.rec=c(}\eqn{\lambda, \alpha}\code{)}.\\
#'   If \code{dist.rec="lognormal"}, the hazard function is
#'   \deqn{\lambda_0(t)=[(1/(\sigma t))*\phi((ln(t)-\mu)/\sigma)]/[\Phi((-ln(t)-\mu)/\sigma)],}
#'   where \eqn{\phi} is the probability density function and \eqn{\Phi} is the cumulative
#'   distribution function of the standard normal distribution, \eqn{\mu\in(-\infty,+\infty)} is a
#'   location parameter and \eqn{\sigma>0} is a shape parameter. Then \code{par.rec=c(}\eqn{\mu,\sigma}\code{)}.
#'   Please note, that specifying \code{dist.rec="lognormal"} together with some covariates does not
#'   specify the usual lognormal model (with covariates specified as effects on the parameters of the
#'   lognormal distribution resulting in non-proportional hazards), but only defines the baseline
#'   hazard and incorporates covariate effects using the proportional hazard assumption.\\
#'   If \code{dist.rec="step"} the hazard function is \deqn{\lambda_0(t)=a, t<=t_1, and \lambda_0(t)=b, t>t_1}.
#'   Then \code{par.rec=c(}\eqn{a,b,t_1}\code{)}.
#' @param pfree Probability that after experiencing an event the individual is not at risk
#'   for experiencing further events for a length of \code{dfree} time units.
#'   Default is \code{pfree=0}.
#' @param dfree Length of the risk-free interval. Must be in the same time unit as \code{fu.max}.
#'  Default is \code{dfree=0}, i.e. the individual is continously at risk for experiencing
#'  events until end of follow-up.

# SIMRECEV - SIMULAÇÃO DE EVENTOS RECORRENTES ====
spsimrec <- function(N,
                      fu.min,
                      fu.max,
                      cens.prob = 0,
                      dist.x = NULL,
                      par.x = NULL,
                      beta.x = NULL,
                      dist.z = c("gamma","lognormal"),
                      random.ef=0,
                      ent.dist.z=0,
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
  nr.cov <- length(beta.x)

  ## gen_fu - Gera períodos de acompanhamento  ====
  gen_fu<-function(N,cens.prob,fu.min,fu.max){
    fu <- rbinom(N, 1, cens.prob) # 1 = censored
    nr.cens <- sum(fu)
    if (nr.cens == 0) { # nobody censored
      fu <- runif(N, min = fu.min, max = fu.max)}
    else {
      index.cens <- which(fu == 1)
      fu[-index.cens] <- runif((N - nr.cens), min = fu.min, max = fu.max)
      fu[index.cens] <- runif(nr.cens, min = 0, max = fu.max)
    }
    return(fu)
  }
  set.seed(123)
  fu<-gen_fu(N,cens.prob,fu.min,fu.max)
  set.seed(NULL)

  ## gen_cov - Gera covariáveis  ====
  gen_cov<-function(N,
                    dist.x,
                    par.x,
                    beta.x){
    nr.cov <- length(beta.x)
    x <- matrix(0, N, nr.cov)
    for (i in 1:nr.cov) {
      dist.x[i] <- match.arg(dist.x[i], choices = c("binomial", "normal"))
      if (dist.x[i] == "binomial") {
        x[, i] <- c(rbinom(N, 1, par.x[[i]]))

      } else { # normally distributed covariate
        mu.x <- par.x[[i]][1]
        sigma.x <- par.x[[i]][2]
        x[, i] <- c(rnorm(N, mean = mu.x, sd = sigma.x))
      }
    }
    return(x)
  }

  set.seed(123)
  if(nr.cov!=0){
    x<-gen_cov(N,
               dist.x,
               par.x,
               beta.x)
    x1<-as.data.frame(cbind(ID,x))
    colnames(x1)<-c("ID","X1","X2")
  }else{x<-rep(0,N)
  x1<-as.data.frame(cbind(ID,x))
  colnames(x1)<-c("ID","X")}
  set.seed(NULL)

  ## gen_rnd_ef - Gera efeitos aleatórios  ====
  gen_rnd_ef<-function(N, ID, dist.z, ent.dist.z, par.z,mu.omega,sigma.omega){
    if (ent.dist.z==0){
      if(par.z==0){# if par.z=0 then frailty=1 for all
        z <- rep(1, N)}
      else{
        dist.z <- match.arg(dist.z, choices = c("gamma", "lognormal"))
        if (dist.z == "gamma") { # gamma-frailty  -- #{Y_i(t) * \lambda_0(t)* Z_i *exp(\beta^t X_i)}
          aGamma <- 1 / par.z
          rnd_ef <- rgamma(N, shape = aGamma, scale = 1 / aGamma)}
        else { # lognormal  -- E(Z)=1
          mu <- log(1 / sqrt(par.z + 1))
          sigma <- sqrt(log(par.z + 1))
          rnd_ef <- exp(rnorm(N, mean = mu, sd = sigma))
        }
      }
    }
    else{ # lognormal  #{Y_i(t) * \lambda_0(t)*exp(\beta^t X_i+\omega_i)}
      rnd_ef <- rnorm(N, mean = mu.omega, sd = sigma.omega)
    }
    return(rnd_ef)
  }

  if(random.ef!=0){
    set.seed(321)
    rnd_ef<- gen_rnd_ef(N, ID, dist.z, ent.dist.z, par.z,mu.omega,sigma.omega)
    set.seed(NULL)
  }else{
    rnd_ef<-rep(1,N)
    ent.dist.z<-0
  }

  rnd_ef1<-as.data.frame(cbind(ID,rnd_ef))
  if(ent.dist.z==0){colnames(rnd_ef1)<-c("ID","z")}
  else{colnames(rnd_ef1)<-c("ID","w")}

  ## Define indivíduos recorrentes (INFLAÇÃO DE ZEROS)  ====

  gen_zi<-function(ID,N,pi){
  recurr <- t(rbinom(N, 1, pi))
  }
  set.seed(234)
  recurr<-gen_zi(ID,N,pi)
  set.seed(NULL)

  recurr1<-as.data.frame(t(rbind(ID,recurr)))
  colnames(recurr1)<-c("ID","recurr")

  # if(logist==1){
  #  pi<-1/(1+exp(-(1+x %*% beta.x)))
  # }


  ## gen_data - Gera tempos de ocorrência dos eventos  ====
  gen_data<-function(ID,
                     N,
                     dist.rec,
                     par.rec,
                     fu,
                     x,
                     rnd_ef){

    if (dist.rec == "weibull") { # weibull
      alpha1 <- par.rec[1]
      alpha2 <- par.rec[2]
    }

    ## Cálculo de alpha1_este e exp_eta ====
    # Considera a forma utilizada para introdução de efeitos aleatórios
    if(nr.cov==0){exp_eta=rep(1,N)}
    else if(ent.dist.z==0){#{Y_i(t) * \lambda_0(t)* Z_i *exp(\beta^t X_i)}
      exp_eta <- exp(x %*% beta.x) * rnd_ef
    }else{#{Y_i(t) * \lambda_0(t)*exp(\beta^t X_i+\omega_i)}
      exp_eta <- exp(x %*% beta.x + rnd_ef)
    }
    alpha1_eta <- alpha1*exp_eta


    ## Definição dos tempos de ocorrência dos primeiros eventos ====
    T<-NULL
    T1<-NULL
    IND<-NULL
    for (i in 1:N) {
      t<-NULL
      U <- runif(1)
      if (dist.rec == "weibull") {
        t <- ((-1)*log(U)*(alpha1_eta[i])^(-1))^(1 / alpha2) #Inversa da função acumulada (veja artigo Generating survival times to simulate pag 1717 tabela II)
        ind<-0
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

  tab <-gen_data(ID, N, dist.rec, par.rec, fu, x,rnd_ef)
  return(tab)
}
