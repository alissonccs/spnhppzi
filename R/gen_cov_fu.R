#' @title gencovfu
#' @aliases gencovfu
#' @export
#' @description
#' Generates covariate values and follow-up periods for recurrent event data simulations.
#' Allows for censoring, different covariate distributions, and regression coefficients.
#'
#' @param N Integer. Number of individuals in the dataset.
#' @param fu.min Numeric. Minimum follow-up time for individuals.
#' @param fu.max Numeric. Maximum follow-up time for individuals.
#' @param cens.prob Numeric. Probability of censoring before the maximum follow-up period.
#' @param dist.x Character. Distribution of covariates (`"binomial"` or `"normal"`).
#' @param par.x List. Parameters for the covariate distributions.
#' @param beta.x Numeric vector. Regression coefficients for the covariates.
#'
#' @return A data frame containing generated covariate values and follow-up times.
#'
#' @examples
#' # Example: Generate covariates with a normal distribution and a binomial variable
#' cov_data <- gencovfu(
#'   N = 500,
#'   fu.min = 7,
#'   fu.max = 7,
#'   cens.prob = 0,
#'   dist.x = c("binomial", "normal"),
#'   par.x = list(0.7, c(0, 1)),
#'   beta.x = c(0.5, 1.3)
#' )
#' head(cov_data)
#'
#' @export

# gencovfu - Gera valores das covariáveis e período de acompanhamento ====
gencovfu <- function(N,
                      fu.min,
                      fu.max,
                      cens.prob = 0,
                      dist.x = NULL,
                      par.x = NULL,
                      beta.x = NULL
                      # int.logist = "FALSE"
                      ){
  ID <- c(1:N)
  # dist.z <- tolower(dist.z)
  # dist.z <- match.arg(dist.z)
  nr.cov <- length(beta.x)

  ##  Gera períodos de acompanhamento  ====
  ### Funçõa - gen_fu ====
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
    fu<-as.data.frame(cbind(ID,fu))
    colnames(fu)<-c("ID","fu")
    return(fu)
  }
 ### Executa função ====
  set.seed(123)
  fu<-gen_fu(N,cens.prob,fu.min,fu.max)
  set.seed(NULL)


  # Gera covariáveis  ====
  ### Função - gen_cov ====
  gen_cov<-function(N,
                    dist.x,
                    par.x,
                    beta.x){
    nr.cov <- length(beta.x)
    x <- matrix(0, N, nr.cov)
    for (i in 1:nr.cov) {
      dist.x[i] <- match.arg(dist.x[i], choices = c("binomial", "normal","intercept.log"))
      if (dist.x[i] == "binomial") {
        x[, i] <- c(rbinom(N, 1, par.x[[i]]))

      } else {
        if (dist.x[i] == "normal") {
        mu.x <- par.x[[i]][1]
        sigma.x <- par.x[[i]][2]
        x[, i] <- c(rnorm(N, mean = mu.x, sd = sigma.x))
        }
        else{
          x[, i] <- rep(1,N)
        }
      }
    }
    return(x)
  }

### Executa Função ====
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

