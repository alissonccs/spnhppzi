# Generated by rstantools.  Do not edit by hand.

# names of stan models
stanmodels <- c("NHPP", "NHPP_RE", "NHPP_RE_1", "NHPP_SE", "SZI_NHPP", "SZI_NHPP_RE", "SZI_NHPP_RE_COV", "SZI_NHPP_SE", "SZI_NHPP_SE_COV", "S_NHPP", "S_NHPP_RE", "S_NHPP_SE", "ZI_NHPP", "ZI_NHPP_RE", "ZI_NHPP_RE_COV", "ZI_NHPP_SE", "ZI_NHPP_SE_COV")

# load each stan module
Rcpp::loadModule("stan_fit4NHPP_mod", what = TRUE)
Rcpp::loadModule("stan_fit4NHPP_RE_mod", what = TRUE)
Rcpp::loadModule("stan_fit4NHPP_RE_1_mod", what = TRUE)
Rcpp::loadModule("stan_fit4NHPP_SE_mod", what = TRUE)
Rcpp::loadModule("stan_fit4SZI_NHPP_mod", what = TRUE)
Rcpp::loadModule("stan_fit4SZI_NHPP_RE_mod", what = TRUE)
Rcpp::loadModule("stan_fit4SZI_NHPP_RE_COV_mod", what = TRUE)
Rcpp::loadModule("stan_fit4SZI_NHPP_SE_mod", what = TRUE)
Rcpp::loadModule("stan_fit4SZI_NHPP_SE_COV_mod", what = TRUE)
Rcpp::loadModule("stan_fit4S_NHPP_mod", what = TRUE)
Rcpp::loadModule("stan_fit4S_NHPP_RE_mod", what = TRUE)
Rcpp::loadModule("stan_fit4S_NHPP_SE_mod", what = TRUE)
Rcpp::loadModule("stan_fit4ZI_NHPP_mod", what = TRUE)
Rcpp::loadModule("stan_fit4ZI_NHPP_RE_mod", what = TRUE)
Rcpp::loadModule("stan_fit4ZI_NHPP_RE_COV_mod", what = TRUE)
Rcpp::loadModule("stan_fit4ZI_NHPP_SE_mod", what = TRUE)
Rcpp::loadModule("stan_fit4ZI_NHPP_SE_COV_mod", what = TRUE)

# instantiate each stanmodel object
stanmodels <- sapply(stanmodels, function(model_name) {
  # create C++ code for stan model
  stan_file <- if(dir.exists("stan")) "stan" else file.path("inst", "stan")
  stan_file <- file.path(stan_file, paste0(model_name, ".stan"))
  stanfit <- rstan::stanc_builder(stan_file,
                                  allow_undefined = TRUE,
                                  obfuscate_model_name = FALSE)
  stanfit$model_cpp <- list(model_cppname = stanfit$model_name,
                            model_cppcode = stanfit$cppcode)
  # create stanmodel object
  methods::new(Class = "stanmodel",
               model_name = stanfit$model_name,
               model_code = stanfit$model_code,
               model_cpp = stanfit$model_cpp,
               mk_cppmodule = function(x) get(paste0("rstantools_model_", model_name)))
})
