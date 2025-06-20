##########################################################################################################################
#' @title spnhppzi: Hierarchical Models for Recurrent Event Data with Zero Inflation and Spatial Correlation
#'
#' @description
#' The `fit_spnhppzi` function estimates parameters for **recurrent event data with zero inflation and spatial correlation**,
#' using either a frequentist or Bayesian approach.
#'
#' @param formula An object of class `"formula"` (or coercible to one).
#'   Defines the symbolic representation of the model to be fitted.
#' @param data A data frame containing all model variables, with the following elements:
#'   - `id`: A vector identifying individuals.
#'   - `status`: A vector indicating whether the event occurred (`1`) or was censored (`0`).
#'   - `stop`: A numeric vector representing the time until the event occurs.
#'   - `IndRec`: A binary indicator (`0/1`) for whether the individual has at least one recurrence.
#'   - `covar`: A vector or matrix of covariates.
#' @param baseline Character. Specifies the baseline intensity function:
#'   - `"plp"`: Parametric Power-Law Process (PLP).
#'   - `"bp"`: Semiparametric approach using Bernstein polynomials.
#' @param approach An integer specifying the estimation approach:
#'   - `mle`: Frequentist approach (available only for models **without** random effects).
#'   - `bayes`: Bayesian approach.
#' @param n_iter An integer specifying the number of iterations for the MCMC algorithm. Default is `4000`.
#' @param n_cores An integer indicating the number of CPU cores to use for parallel computation. Default is `1` (no parallelization).
#' @param n_chains An integer defining the number of Markov chains to run in parallel. Default is `4`.
#' @param ZI A character string indicating whether the model accounts for zero-inflated data:
#'   - `"true"`: The model includes a structure for zero inflation.
#'   - `"false"`: The model does not consider zero inflation.
#' @param rnd_efc  A character string indicating whether the model includes a random effects structure (frailty term):
#'   - `"true"`: Includes random effects (with frailty).
#'   - `"false"`: No random effects (no frailty).
#' @param initial A vector of initial parameter values for the optimization procedure.
#' @param shp_alpha1 Shape parameter of the Gamma prior distribution for the scale parameter (`alpha_1`)
#'   when using the Power-Law Process (PLP) model. Default is `0.1`.
#' @param scl_alpha1 Scale parameter of the Gamma prior distribution for the scale parameter (`alpha_1`)
#'   when using the PLP model. Default is `0.1`.
#' @param shp_alpha2 Shape parameter of the Gamma prior distribution for the shape parameter (`alpha_2`)
#'   when using the PLP model. Default is `0.1`.
#' @param scl_alpha2 Scale parameter of the Gamma prior distribution for the shape parameter (`alpha_2`)
#'   when using the PLP model. Default is `0.1`.
#' @param mu_beta Mean (`μ`) of the Normal prior distribution for the regression coefficients of the covariates. Default is `0`.
#' @param sigma_beta Standard deviation (`σ`) of the Normal prior distribution for the regression coefficients of the covariates. Default is `4`.
#' @param mu_omega Mean (`μ_omega`) of the Normal distribution for individual-level random effects (`omega`)
#'   in models without spatial random effects. Default is `0`.
#' @param shp_sigma2_z Shape parameter of the Gamma prior distribution for the hyperparameter `sigma2_z`,
#'   which controls the variance of the random effects. Default is `0.01`.
#' @param scl_sigma2_z Scale parameter of the Gamma prior distribution for the hyperparameter `sigma2_z`.
#'   Default is `0.01`.
#' @param spatial A logical value indicating whether the model includes a spatial random effects structure:
#'   - `0`: No spatial random effects.
#'   - `1`: Includes spatial random effects.#'
#'   If `spatial = TRUE`, then `rnd_efc` must necessarily be `FALSE`.
#' @param sp_model A character string specifying the spatial model to be used (only applicable if `spatial = TRUE`):
#'   - `"icar"`: Intrinsic Conditional Autoregressive (ICAR) model.
#'   - `"car"`: Conditional Autoregressive (CAR) model (under development).
#'   - `"sparse"`: Sparse spatial model (under development).
#' @param nb_mat A neighborhood matrix defining the spatial structure of the area units.
#'   It is an object of class `sf`, typically representing adjacency relationships between spatial units.
#' @param shp_tau Shape parameter of the Gamma prior distribution for the precision (`tau`) of the spatial random effects.
#'   Default is `0.01`.
#' @param scl_tau Scale parameter of the Gamma prior distribution for the precision (`tau`) of the spatial random effects.
#'   Default is `0.01`.
#' @param W_n An integer specifying the number of edges (connections) in the sparse neighborhood matrix.
#'   It represents the number of area pairs that share a spatial connection, as defined in `W_sparse`.
#'   This value can be computed as:
#'   - If `W` is a binary adjacency matrix: `W_n = sum(W) / 2`.
#'   - If `nb` is a neighborhood list (e.g., from `spdep::poly2nb()`): `W_n = sum(sapply(nb, length)) / 2`.
#'   - If `W_sparse` is a list of connected area pairs: `W_n = nrow(W_sparse)`.
#' @param bp_degree An integer specifying the degree of the Bernstein polynomials used in the semiparametric
#'   baseline intensity function. Default is `NULL`, meaning it is automatically selected.
#'   This argument is only used when `baseline = "bp"`.
#' @param h1_gamma A numeric value specifying the mean of the lognormal prior distribution
#'   for the coefficients of the Bernstein polynomials. Default is `0`.
#'   This argument is only used when `baseline = "bp"`.
#' @param h2_gamma A numeric value specifying the variance of the lognormal prior distribution
#'   for the coefficients of the Bernstein polynomials. Default is `4`.
#'   This argument is only used when `baseline = "bp"`.
#' @return A list containing the estimated parameters and model fit statistics.
#'
#' @export
#'
#' @examples
#' # EXAMPLE 01 ----
#' # This example illustrates the simplest case: the NHPP model for recurrent event data
#'
#' N <- 500
#' alpha1_r <- 0.5
#' alpha2_r <- 1.3
#' beta1_r <- 0.6
#' beta2_r <- 0.8
#' pi_r <- 0
#' fu.min <- 7
#' fu.max <- 7
#'
#' set.seed(5832)
#' cov.fu <- gencovfu(N = N,
#'                     fu.min = fu.min,
#'                     fu.max = fu.max,
#'                     cens.prob = 0,
#'                     dist.x = c("binomial", "normal"),
#'                     par.x = list(0.7, c(0, 1)),
#'                     beta.x = c(beta1_r, beta2_r))
#' set.seed(NULL)
#'
#' base <- spsimrec(N = cov.fu$N,
#'                  cov_rec = c("ID", "X1", "X2"),
#'                  beta_x_rec = c(beta1_r, beta2_r),
#'                  logist = 0,
#'                  x1 = cov.fu$x1,
#'                  fu = cov.fu$fu,
#'                  fu_max = cov.fu$fu.max,
#'                  fu_min = cov.fu$fu.min,
#'                  spatial = 0,
#'                  random_ef = 0,
#'                  pi = pi_r,
#'                  par_z = 0,
#'                  dist_int_func = "weibull",
#'                  par_int_func = c(alpha1_r, alpha2_r),
#'                  baseline = "plp")
#'
#' formula2 <- as.list(Formula(spnhppzi::Recur(time = end, event = status, id = ID, SP_ID = NULL, IndRec = IndRec) ~ X1 + X2 | -1))
#'
#' RESULT_BAYES_SCOV1 <- spnhppzi::fit_spnhppzi(formula2,
#'                                          base,
#'                                          baseline = "plp",
#'                                          rnd_efc = FALSE,
#'                                          ZI = FALSE,
#'                                          approach = "BAYES",
#'                                          sp_model = "ICAR",
#'                                          initial = 1,
#'                                          shp_alpha1 = 0.1, scl_alpha1 = 0.1,
#'                                          shp_alpha2 = 0.1, scl_alpha2 = 0.1,
#'                                          mu_beta = 0, sigma_beta = 4,
#'                                          mu_psi = 0, sigma_psi = 4,
#'                                          mu_omega = 0,
#'                                          spatial = 0,
#'                                          n_iter = 2000,
#'                                          n_cores = 2,
#'                                          n_chains = 2)
#'
#' summary(RESULT_BAYES_SCOV1, pars = c("alpha", "beta"))
#'
#' # EXAMPLE 02 ----
#' # This example illustrates the SZI-NHPP-SE model (Semiparametric Zero-Inflated Non-Homogeneous Poisson Process with Spatial Effects)
#'
#' # Load adjacency matrix from package's extdata directory
#' Adj_matrix <- readRDS(system.file("extdata", "Adj_matrix.RDS", package = "spnhppzi"))
#' list_area_RMBH<-as.numeric(row.names(Adj_matrix))
#' # Define parameters
#' N <- 500
#' alpha1_r <- 2
#' alpha2_r <- 1.3
#' beta1_r <- 0.6
#' beta2_r <- 0.8
#' sp_tau_r <- 1
#' psi1_r <- 1.6
#' psi2_r <- 1.2
#' pi_r <- 0.75
#' fu.min <- 7
#' fu.max <- 7
#' degree_bp <- min(ceiling(N^0.4), 5)
#'
#' # Simulating covariates for the dataset
#' cov.fu <- gencovfu(
#'   N = N,
#'   fu.min = fu.min,
#'   fu.max = fu.max,
#'   cens.prob = 0,
#'   dist.x = c("binomial", "normal"),
#'   par.x = list(0.7, c(0, 1)),
#'   beta.x = c(beta1_r, beta2_r)
#' )
#'
#' # Simulating recurrent event data for model estimation
#' base_sp <- spsimrec(
#'   N = cov.fu$N,
#'   cov_rec = c("ID", "X1", "X2"),
#'   beta_x_rec = c(beta1_r, beta2_r),
#'   logist = 0,
#'   x1 = cov.fu$x1,
#'   fu = cov.fu$fu,
#'   fu_max = cov.fu$fu.max,
#'   fu_min = cov.fu$fu.min,
#'   spatial = 1,
#'   list_area = list_area_RMBH,
#'   sp_model = "ICAR",
#'   SP_N = 133,
#'   nb_mat = Adj_matrix,
#'   sp_tau = sp_tau_r,
#'   random_ef = 1,
#'   pi = pi_r,
#'   par_z = 0,
#'   dist_int_func = "weibull",
#'   par_int_func = c(alpha1_r, alpha2_r),
#'   baseline = "plp"
#' )
#'
#' # Fitting the SZI-NHPP-SE model
#' formula2 <- as.list(Formula(spnhppzi::Recur(end, status, ID, SP_ID, IndRec) ~ X1 + X2 | -1))
#'
#' RESULT <- spnhppzi::fit_spnhppzi(
#'   formula2,
#'   base_sp,
#'   baseline = "bp",
#'   rnd_efc = TRUE,
#'   ZI = TRUE,
#'   approach = "BAYES",
#'   sp_model = "ICAR",
#'   initial = 1,
#'   mu_beta = 0, sigma_beta = 4,
#'   mu_psi = 0, sigma_psi = 4,
#'   mu_omega = 0,
#'   spatial = 1,
#'   nb_mat = Adj_matrix,
#'   shp_tau = 0.01,
#'   scl_tau = 0.01,
#'   n_iter = 2000,
#'   n_cores = 1,
#'   n_chains = 2,
#'   W_n = 365,
#'   bp_degree = degree_bp,
#'   h1_gamma = 0,
#'   h2_gamma = 4
#' )
#'
#' summary(RESULT$result_stan, pars = c("alpha", "beta", "pii", "tau"))#'
#'
#' #Example 03: ZI-NHPP-SE Model (Real Data Application) ----
#' #This example applies the Zero-Inflated Non-Homogeneous Poisson Process with Spatial Effects (SZI-NHPP-SE) model to real criminal recidivism data, as presented in the application section of the paper accepted in JRSS-A.
#' #The model includes spatially structured random effects and handles excess zeros. The covariate sex was used to assess differences in recidivism patterns.
#'
#' # Load the bodily_injury dataset from the package's extdata directory
#' df_bodily_injury <- readRDS(system.file("extdata", "df_bodily_injury.RDS", package = "spnhppzi"))
#'
#' # Load adjacency matrix from package's extdata directory
#' Adj_matrix_aplication <- readRDS(system.file("extdata", "Adj_matrix_aplication.RDS", package = "spnhppzi"))
#'
#' # Fitting the SZI-NHPP-SE model
#' formula2=as.list(Formula(spnhppzi::Recur(time=end,event=as.numeric(status),id=id1,SP_ID=SP_ID,IndRec=IndRec)~sexo1|-1))
#' RESULT_BAYES_SCOV1<- spnhppzi::fit_spnhppzi(formula2,
#'                                            df_bodily_injury,
#'                                            baseline = "plp",
#'                                            rnd_efc = "TRUE",
#'                                            ZI="TRUE",
#'                                            approach = "BAYES",
#'                                            sp_model = "ICAR",
#'                                            initial=1,
#'                                            shp_alpha1=0.1,scl_alpha1=0.1,
#'                                            shp_alpha2=0.1,scl_alpha2=0.1,
#'                                            mu_beta=0,sigma_beta=4,
#'                                            mu_psi=0,sigma_psi=4,
#'                                            mu_omega=0,
#'                                            spatial=1,
#'                                            nb_mat=Adj_matrix_aplication,
#'                                            shp_tau=0.01,
#'                                            scl_tau=0.01,
#'                                            n_iter = 2000,
#'                                            n_cores=1,
#'                                            n_chains=2,
#'                                            W_n=365
#' )
#' RESULT_BAYES_SCOV1_sp_plp<-RESULT_BAYES_SCOV1
#' pars_desc_sp_plp<-summary(RESULT_BAYES_SCOV1_sp_plp, pars = c("alpha","beta","pii","tau"))
#' pars_desc_sp_plp$summary
##########################################################################################################################
fit_spnhppzi<-function(formula,
                       data,
                       baseline = c("plp","bp"),
                       approach = c("mle", "bayes"),
                       n_iter=4000,
                       n_cores=1,
                       n_chains=4,
                       ZI = c("true","false"),
                       rnd_efc= c("true","false"),
                       initial,
                       shp_alpha1=0.1,scl_alpha1=0.1,shp_alpha2=0.1,scl_alpha2=0.1,
                       mu_beta=0,sigma_beta=4,
                       mu_psi=0,sigma_psi=4,
                       mu_omega=0,
                       shp_sigma2_z=0.01, scl_sigma2_z=0.01,
                       spatial=0,
                       sp_model = c("car","sparse","icar"),
                       nb_mat=NULL,
                       shp_tau=0.01,
                       scl_tau=0.01,
                       W_n=0,
                       bp_degree=NULL,
                       h1_gamma=0,
                       h2_gamma=4
                       ){

  formula <- Formula::Formula(formula)
  approach <- tolower(approach)
  approach <- match.arg(approach)
  baseline <- tolower(baseline)
  baseline <- match.arg(baseline)
  ZI <- tolower(ZI)
  ZI <- match.arg(ZI)
  rnd_efc <- tolower(rnd_efc)
  rnd_efc <- match.arg(rnd_efc)
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
                     "plp" = 1,
                     "bp" = 2
  )
  if(baseline == 1){
    m <- 2
  }

  g<-NULL
  G<-NULL
  if(baseline == 2){
    bases <- bp(time, max_stop, bp_degree, zeta, N, n)
    g <- bases$b
    G <- bases$B
    m <- bp_degree
  }

  approach <- switch(approach,
                     "mle" = 0,
                     "bayes" = 1
  )

  rnd_efc <- switch(rnd_efc,
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
                     approach=approach, rnd_efc=rnd_efc, ZI=ZI, begin_ind=begin_ind,end_ind=end_ind,
                     n_ind=n_ind,n_ind1=n_ind1, m=m, mu_omega=mu_omega, shp_sigma2_z=shp_sigma2_z,
                     scl_sigma2_z=scl_sigma2_z, shp_alpha1=shp_alpha1, scl_alpha1=scl_alpha1,
                     shp_alpha2=shp_alpha2, scl_alpha2=scl_alpha2,mu_beta=mu_beta, sigma_beta=sigma_beta,
                     mu_psi=mu_psi, sigma_psi=sigma_psi, baseline=baseline, SP_ID,gr_SP_ID, SP_N, nb_mat=nb_mat,
                     shp_tau=shp_tau, scl_tau=scl_tau, W_n=W_n, G=G, g=g, zeta=zeta
                     )
  # Modelos não espaciais ----
  if(spatial==0){
  ## Sem efeito aleatório para o indivíduo ----
  if(rnd_efc==0){
    ### Semiparamétrico ----
    if(baseline==2){
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
      if(baseline==2){
       mod <-stanmodels$S_NHPP_RE
      }
      #### Paramétrico ----
      else{
      # if(tp_DIC==0){#Usa Transformed Parameters
      # mod <- stanmodels$NHPP_RE
      # }
      # else{#Usa Generate Quantities
      mod <- stanmodels$NHPP_RE_1
      # }
     }
    }
    #### Com inflaćão de zeros ----
    else{
      # Sem covariável na reg logística ---
      if(q==0){
          # Semiparamétrico
          if(baseline==2){
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
        if(baseline==2){
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
    if(baseline==2){
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
          if(baseline==2){
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
        if(baseline==2){
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
  if(rnd_efc==0){
    #Optimizing
    if(approach==0){
      result<- optimizing(stanmodels$mod, data = data_model, hessian =TRUE, init=initial, algorithm = "LBFGS")
      if(baseline==2){
       result<-list("result_stan"=result,"G"=G, "g"=g, "m"=m)
      }
      return(result)
    }
    if(approach==1){
     result_b<- sampling(mod, data = data_model, cores = n_cores, iter=n_iter, chains=n_chains)
      if(baseline==2){
       result_b<-list("result_stan"=result_b,"G"=G, "g"=g, "m"=m)
      }
      return(result_b)
    }
  }
  else{
   result_c<- sampling(mod, data = data_model, cores = n_cores, iter=n_iter,chains=n_chains)
   if(baseline==2){
    result_c<-list("result_stan"=result_c,"G"=G, "g"=g, "m"=m)
   }
    return(result_c)
  }
}
