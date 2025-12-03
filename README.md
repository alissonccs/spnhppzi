# Installation

You can install The spnhppzi package from GitHub using devtools:

```r
## Install devtools if not already installed
install.packages("devtools")

## Install spnhppzi from GitHub
devtools::install_github("alissonccs/spnhppzi@master")
```
After installation, load the package using:

```r
library(spnhppzi)
```

# **spnhppzi: Bayesian Modeling of Recurrent Event Data with Zero Inflation and Spatial Correlation**

The **spnhppzi** package (**Spatial Non-Homogeneous Poisson Process with Zero Inflation**) is designed for modeling recurrent event data that exhibit **zero inflation and spatial correlation**. It follows a **Bayesian approach**, implementing **hierarchical models** to capture **complex structures** associated with event repetition and spatial dependencies.

Spatial modeling is based on the **Intrinsic Conditional Autoregressive (ICAR) model**, which efficiently incorporates spatial correlation. Additionally, the package provides **both parametric and semiparametric models**, where **Bernstein polynomials** are used for modeling the **baseline intensity function**. This approach enhances the model's flexibility, making it applicable in scenarios where **purely parametric models may struggle to adequately capture data complexity**.

## **Data Simulation**

The **spnhppzi** package also includes **functions for simulating recurrent event data**, using an adaptation of the **SIMREC** package ([Farrington et al., 2014](https://cran.r-project.org/web/packages/simrec/vignettes/simrec-vignette.html)). This enables the generation of **spatially correlated recurrent event data**, allowing for **model performance evaluation** and experimentation with different recurrence and spatial dependence scenarios.

## **References and Application**

This repository contains the functions used to generate the results presented in the following article:

📄 **"The Analysis of Criminal Recidivism: A Hierarchical Model-Based Approach for the Analysis of Zero-Inflated, Spatially Correlated Recurrent Events Data"**, available at [Journal of the Royal Statistical Society - Series A]([https://arxiv.org/abs/2405.02666](https://academic.oup.com/jrsssa/advance-article/doi/10.1093/jrsssa/qnaf061/8152038)).

Additional methodological details can be found in the dissertation:

📖 **"Hierarchical Models for the Analysis of Recurrent Event Data with Zero Inflation and Spatial Correlation"**, available upon request from the author via **alisson.ccs2@gmail.com**.

## **Features**

- **Bayesian hierarchical modeling** for recurrent event data;
- **Spatial correlation modeling** using the **ICAR model**;
- **Parametric and semiparametric models**, incorporating **Bernstein polynomials** for greater flexibility in modeling the baseline intensity function;
- **Data simulation** for spatially correlated recurrent events (adaptation of **SIMREC**);
- **Tools for data preprocessing and model estimation**.

## **Examples**

### **Example 1: NHPP Model for Recurrent Event Data**
```r
# EXAMPLE ----
# This example illustrates the simplest case: the NHPP model for recurrent event data

N <- 500
alpha1_r <- 0.5
alpha2_r <- 1.3
beta1_r <- 0.6
beta2_r <- 0.8
pi_r <- 0
fu.min <- 7
fu.max <- 7

set.seed(5832)
cov.fu <- gencovfu(N = N,
                     fu.min = fu.min,
                     fu.max = fu.max,
                     cens.prob = 0,
                     dist.x = c("binomial", "normal"),
                     par.x = list(0.7, c(0, 1)),
                     beta.x = c(beta1_r, beta2_r))
set.seed(NULL)

base <- spsimrec(N = cov.fu$N,
                  cov_rec = c("ID", "X1", "X2"),
                  beta_x_rec = c(beta1_r, beta2_r),
                  logist = 0,
                  x1 = cov.fu$x1,
                  fu = cov.fu$fu,
                  fu_max = cov.fu$fu.max,
                  fu_min = cov.fu$fu.min,
                  spatial = 0,
                  random_ef = 0,
                  pi = pi_r,
                  par_z = 0,
                  dist_int_func = "weibull",
                  par_int_func = c(alpha1_r, alpha2_r),
                  baseline = "plp")

formula2 <- as.list(Formula(spnhppzi::Recur(time = end, event = status, id = ID, SP_ID = NULL, IndRec = IndRec) ~ X1 + X2 | -1))

RESULT_BAYES_SCOV1 <- spnhppzi::fit_spnhppzi(formula2,
                                          base,
                                          baseline = "plp",
                                          rnd_efc = FALSE,
                                          ZI = FALSE,
                                          approach = "BAYES",
                                          sp_model = "ICAR",
                                          initial = 1,
                                          shp_alpha1 = 0.1, scl_alpha1 = 0.1,
                                          shp_alpha2 = 0.1, scl_alpha2 = 0.1,
                                          mu_beta = 0, sigma_beta = 4,
                                          mu_psi = 0, sigma_psi = 4,
                                          mu_omega = 0,
                                          spatial = 0,
                                          n_iter = 2000,
                                          n_cores = 2,
                                          n_chains = 2)

summary(RESULT_BAYES_SCOV1, pars = c("alpha", "beta"))
```

### **Example 2: SZINHPP Model with Spatial Correlation**
```r
# ADDITIONAL EXAMPLE ----
# This example illustrates the SZINHPP model (Spatial Zero-Inflated NHPP) with spatial correlation

# Load adjacency matrix from package's extdata directory
Adj_matrix <- readRDS(system.file("extdata", "Adj_matrix.RDS", package = "spnhppzi"))
list_area_RMBH <- as.numeric(row.names(Adj_matrix))

# Define parameters
N <- 500
alpha1_r <- 2
alpha2_r <- 1.3
beta1_r <- 0.6
beta2_r <- 0.8
sp_tau_r <- 1
psi1_r <- 1.6
psi2_r <- 1.2
pi_r <- 0.75
fu.min <- 7
fu.max <- 7
degree_bp <- min(ceiling(N^0.4), 5)

# Simulating covariates for the dataset
cov.fu <- gencovfu(
  N = N,
  fu.min = fu.min,
  fu.max = fu.max,
  cens.prob = 0,
  dist.x = c("binomial", "normal"),
  par.x = list(0.7, c(0, 1)),
  beta.x = c(beta1_r, beta2_r)
)

# Simulating recurrent event data for model estimation
base_sp <- spsimrec(
  N = cov.fu$N,
  cov_rec = c("ID", "X1", "X2"),
  beta_x_rec = c(beta1_r, beta2_r),
  logist = 0,
  x1 = cov.fu$x1,
  fu = cov.fu$fu,
  fu_max = cov.fu$fu.max,
  fu_min = cov.fu$fu.min,
  spatial = 1,
  list_area = list_area_RMBH,
  sp_model = "ICAR",
  SP_N = 133,
  nb_mat = Adj_matrix,
  sp_tau = sp_tau_r,
  random_ef = 1,
  pi = pi_r,
  par_z = 0,
  dist_int_func = "weibull",
  par_int_func = c(alpha1_r, alpha2_r),
  baseline = "plp"
)
# Fitting the SZINHPP model
formula2 <- Formula(spnhppzi::Recur(end, status, ID, SP_ID, IndRec) ~ X1 + X2 | -1)
RESULT <- spnhppzi::fit_spnhppzi(
  formula2,
  base_sp,
  baseline = "bp",
  rnd_efc = TRUE,
  ZI = TRUE,
  approach = "BAYES",
  sp_model = "ICAR",
  initial = 1,
  mu_beta = 0, sigma_beta = 4,
  mu_psi = 0, sigma_psi = 4,
  mu_omega = 0,
  spatial = 1,
  nb_mat = Adj_matrix,
  shp_tau = 0.01,
  scl_tau = 0.01,
  n_iter = 2000,
  n_cores = 1,
  n_chains = 2,
  W_n = 365,
  bp_degree = degree_bp,
  h1_gamma = 0,
  h2_gamma = 4
)

summary(RESULT$result_stan, pars = c("alpha", "beta", "pii", "tau"))
```
### **Example 3: ZI-NHPP-SE Model for Real Data Application**

This example illustrates the application of the **Zero-Inflated Non-Homogeneous Poisson Process with Spatial Effects (ZI-NHPP-SE)** model to real criminal recidivism data.

The model structure and its implementation are based on the methodology presented in the paper:  
📄 *“The Analysis of Criminal Recidivism: A Hierarchical Model-Based Approach for the Analysis of Zero-Inflated, Spatially Correlated Recurring Event Data”*, accepted for publication in the **Journal of the Royal Statistical Society: Series A (Statistics in Society)**.

In this example, we consider **criminal recidivism events** across different regions, accounting for **spatial dependence** among areas and **excess zeros** (individuals with no recidivism). The model also includes the individual-level covariate sex, used to assess differences in recidivism between males and females.
The code below reproduces the **ZI-NHPP-SE** model used in the empirical application section of the article.

```r
# Load the bodily_injury dataset from the package's extdata directory
df_bodily_injury <- readRDS(system.file("extdata", "df_bodily_injury.RDS", package = "spnhppzi"))

# Load adjacency matrix from package's extdata directory
Adj_matrix_aplication <- readRDS(system.file("extdata", "Adj_matrix._aplication.RDS", package = "spnhppzi"))

# Build the formula object
formula2 <- as.list(Formula(
  spnhppzi::Recur(time = end,
                  event = as.numeric(status),
                  id = id1,
                  SP_ID = SP_ID,
                  IndRec = IndRec) ~ sexo1 | -1))

# Fit the model
RESULT_BAYES_SCOV1 <- spnhppzi::fit_spnhppzi(
  formula2,
  df_bodily_injury,
  baseline = "plp",
  rnd_efc = TRUE,
  ZI = TRUE,
  approach = "BAYES",
  sp_model = "ICAR",
  initial = 1,
  shp_alpha1 = 0.1, scl_alpha1 = 0.1,
  shp_alpha2 = 0.1, scl_alpha2 = 0.1,
  mu_beta = 0, sigma_beta = 4,
  mu_psi = 0, sigma_psi = 4,
  mu_omega = 0,
  spatial = 1,
  nb_mat = Adj_matrix_aplication,
  shp_tau = 0.01,
  scl_tau = 0.01,
  n_iter = 2000,
  n_cores = 1,
  n_chains = 2,
  W_n = 365
)

# Extract posterior summaries
RESULT_BAYES_SCOV1_sp_plp <- RESULT_BAYES_SCOV1
pars_desc_sp_plp <- summary(RESULT_BAYES_SCOV1_sp_plp, pars = c("alpha", "beta", "pii", "tau", "omega"))
pars_desc_sp_plp <- pars_desc_sp_plp$summary



