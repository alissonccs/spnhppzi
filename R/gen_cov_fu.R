#'@title gencovfu
#' @aliases gencovfu
#' @export
#' @description Recur:
#' @param cens.prob  Probabiidade de apresentar censura antes do período máximo de acompanhamento.
#' @param dist.x     Distribuição das Covariáveis. Binomial ou normal.
#' @param par.x      Parâmetros das distribuições das covariáveis.
#' @param beta.x     Coeficiente de regressão das covariáveis.

# gencovfu - Gera valores das covariáveis e período de acompanhamento ====
gencovfu <- function(N,
                      fu.min,
                      fu.max,
                      cens.prob = 0,
                      dist.x = NULL,
                      par.x = NULL,
                      beta.x = NULL
                      ){
  ID <- c(1:N)
  # dist.z <- tolower(dist.z)
  # dist.z <- match.arg(dist.z)
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

      } else {
        mu.x <- par.x[[i]][1]
        sigma.x <- par.x[[i]][2]
        x[, i] <- c(rnorm(N, mean = mu.x, sd = sigma.x))
      }
    }
    return(x)
  }

  set.seed(1)
  if(nr.cov!=0){
    x<-gen_cov(N,
               dist.x,
               par.x,
               beta.x)
    x1<-as.data.frame(cbind(ID,x))
    colnames(x1)<-c("ID",paste("X",c(1:nr.cov),sep=""))
  }else{x<-rep(0,N)
  x1<-as.data.frame(cbind(ID,x))
  colnames(x1)<-c("ID","X")}
  set.seed(NULL)

  #cov.fu<-cbind(x,fu)
  return(list(N=N,x=x,x1=x1,fu=fu,nr.cov=nr.cov,beta.x=beta.x,fu.max=fu.max,fu.min=fu.min))
}

