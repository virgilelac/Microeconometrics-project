set.seed(1234)

options(error = NULL) #Avoids having the (annoying) R debugging tool at each error 
options(rlang_backtrace_on_error = "none")

install.packages(c("MASS", "parallel", "dplyr", "gmm", "Matrix"))
library(MASS)
library(parallel)
library(dplyr)
library(gmm)
library(Matrix)
library(ggplot2)
#### Replication code for the microeconometrics project. 

##Low level functions (for easier readability when you correct)

groupmaker <- function(N,G){
  ###This function creates arbitrary clusters based on the nuber of clusters and the size of the population.
  ###It can handle uneven cluster sizes 
  ##Args 
  
  #Inputs: 
  #N: number of observations
  #G: number of clusters
  
  #Output: 
  #id: a vector of group assignments
  
  n <- N%/%G
  add <- N%%G
  
  grp_sizes <- c(rep(n+1,add), rep(n, G - add))
  
  if(length(grp_sizes) != G){stop("The vector of group assignments is not of the same length as the specified number of groups.")}
  if(sum(grp_sizes) != N){stop("The number of individuals assigned to a group is smaller than the specified number of individuals.")}
  id <- rep(0,N)
  
  id <- rep(seq_len(G), times = grp_sizes)
}

covariates_gen <- function(id, N, p, G,clustsd){
  ###This function generates cross-sectionnaly dependant covariates. 
  ##Args:
  
  #Inputs: 
  #id: vector of cluster memberships
  #N: number of individuals
  #p: number of covariates
  #clustsd: standard deviation of the ccommon characteristics within the cluster
  
  #Output: a N*p matrix of covariates
  
  grp_dummies <- model.matrix(~factor(id) - 1) #NxG
  common_chars <- matrix(rnorm(G* p, sd = clustsd), G, p) #G*p
  
  X <- grp_dummies %*% common_chars
  
  return(X)
}

varcov_ins <- function(l, rho, sigmaZ) {
  ### This function creates adaptative variance-covariance matrices for the instruments
  ##Args 
  
  #Inputs: 
  #l: number of instruments
  #rho: parameter for the covariance
  #sigmaZ: variance of the instruments
  
  #Output: a lxl variance-covariance matrix
  
  vcov <- matrix(rho, l, l); diag(vcov) <- 1
  S <- sigmaZ * vcov
  
  return(S)
  
}

g <- function(data,theta) {
  X <- model.matrix(~ X1 + X2 + X3, data = data)              
  Z <- model.matrix(~ 0 + Z1 + Z2 + Z3 + Z4 + X3, data = data)    
  u <- as.numeric(data$Y - X %*% theta)                       
  Z * u                                                       
}

CCE_Hwang <- function(N, G, data, theta){
  
  ## We compute g (expressed as f in Hwang 2021) for each individual 
  
  g_hat <- g(data, theta)
  
  ## We compute \sum_{i = 1}^{L_g} f_i^g(\theta)

  id <- data[[1]]
  L_bar <- (1/G) * sum(as.integer(table(data$id)))
  sum_i <- rowsum(g_hat, group = id)
  sum_i <- (1/(sqrt(L_bar)))*sum_i
  ## We compute CCE 
  
  CCE <- crossprod(sum_i) / nrow(sum_i)
  
}

#################################################################################################
############################ HIGHER LEVEL FUNCIONS ##############################################
#################################################################################################

dgp_clust <-function( 
    N = 10000,                  # Number of observations
    G = 25,                     # Number of clusters
    p = 3,                      # Number of parameters in the parameter vector
    endog = 2,                  # Number of endogenous regressors (first 'endog' columns of X)
    serialcor = 0,              # (unused here)
    clustersd = 1,              # SD of cluster common component in X
    errors_var = 1,             # Var(idiosyncratic outcome error)
    l = 4,                      # Number of instruments
    rho_Z = 0.3,                # Cov parameter for the idiosyncratic instrument component
    sigmaZ = 1,                 # Variance of the idiosyncratic instrument component
    alpha = 0.5,
    rho_xe = 0.6,               # Corr between cluster shocks v_g (X) and u_g (Y)
    indivsd = 1,                # SD of idiosyncratic X noise
    sigmaZg = 1,                # Var of cluster-level instrument component Z_g
    target_R2 = 0.4,            # Target first-stage R^2 for each endogenous regressor
    sigma_u = 1.0,              # Var(cluster outcome shock u_g)
    sigma_v = 1.0,              # Var(cluster regressor shock v_g)
    seed = NULL,
    beta = NULL                 # true beta; if NULL draw N(0,1)
){
  if (!is.null(seed)) set.seed(seed)
  if (l <= p) stop("Number of instruments should exceed the number of parameters")
  if (l < endog) stop("Number of instruments l must be >= number of endogenous regressors 'endog'.")
  
  ### Instruments 
  S_ins <- varcov_ins(l, rho_Z, sigmaZ)
  
  ### Group structure (unbalanced clusters allowed)
  id <- groupmaker(N, G)
  
  ### Cluster shocks
  
  Sigma_uv <- matrix(c(sigma_v^2, rho_xe * sigma_v * sigma_u,
                       rho_xe * sigma_v * sigma_u, sigma_u^2), 2, 2)
  
  uv_g <- mvrnorm(G, mu = c(0, 0), Sigma = Sigma_uv)
  v_g  <- uv_g[, 1] 
  u_g  <- uv_g[, 2]
  
  ### Covariates
  X <- covariates_gen(id, N, p, G, clustersd)
  X <- X + matrix(rnorm(N * p, sd = indivsd), N, p)
  
  # We add cluster schock to endogenous
  for (j in seq_len(endog)) X[, j] <- X[, j] + v_g[id]
  
  ### Instruments generation
  Zg   <- mvrnorm(G, mu = rep(0, l), Sigma = diag(sigmaZg, l))
  zeta <- mvrnorm(N, mu = rep(0, l), Sigma = S_ins)
  Z    <- zeta + Zg[id, , drop = FALSE]
  
  ### First stage
  
  Pi_raw <- matrix(rnorm(l * endog), l, endog)
  for (j in seq_len(endog)) {
    signal_ratio <- target_R2 / (1 - target_R2)
    var_Xj_noise <- var(X[, j])
    var_Zbj      <- var(as.vector(Z %*% Pi_raw[, j]))
    scale_j      <- sqrt(signal_ratio * var_Xj_noise / (var_Zbj + 1e-12))
    Pi_raw[, j]  <- scale_j * Pi_raw[, j]
    X[, j]       <- X[, j] + as.vector(Z %*% Pi_raw[, j])
  }
  
  ### Outcome
  
  if (is.null(beta)) beta <- rnorm(p)
  eps <- rnorm(N, sd = sqrt(errors_var))
  Y   <- alpha + as.vector(X %*% beta) + u_g[id] + eps
  
  ### Moments at θ0 
  resid0    <- as.vector(Y - X %*% beta)  # at true θ0
  M         <- Z * resid0                 # N x l
  Bmg       <- rowsum(M, group = id)      # G x l
  Omega_hat <- crossprod(Bmg) / G         # l x l
  
  ##Constitution of dataset
  
  colnames(X) <- paste0("X", seq_len(ncol(X)))
  colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
  data <- data.frame(id = id, Y, X, Z)
  
  return(list(
    data = data,
    beta = beta,
    alpha = alpha,
    Pi = Pi_raw,
    rho_xe = rho_xe,
    Bmg = Bmg,
    Omega_hat = Omega_hat,
    cluster_sizes = as.integer(table(id)),
    details = list(
      p = p, q = endog, l = l, target_R2 = target_R2,
      variances = c(sigma_u = sigma_u, sigma_v = sigma_v, sigma_e = errors_var,
                    sigmaZ = sigmaZ, sigmaZg = sigmaZg),
      rho_Z = rho_Z
    )
  ))
}

simulation <- function(M = 500, 
                       N = 1000,                  #Number of observations
                       G = 25,                     #Number of clusters
                       serialcor = 0.7,              #To allow for serial autocorrelation
                       clustersd = 1,              #Variance within each cluster (covariates)
                       errors_var = 1,             #Variance within each cluster (error terms)
                       rho_Z = 0.3,                #Parameter for the covariance of the instruments
                       sigmaZ = 1
                       ,level = 0.05){
  
  theta <- rep(0,M)
  bias_intercept <<- rep(0,M)
  bias_b1 <- rep(0,M)
  bias_b2 <- rep(0,M)
  coverage_int <- rep(0,M)
  coverage_b1 <- rep(0,M)
  coverage_b2 <- rep(0,M)
  z <- qnorm(1-level/2)
  
  true_b1 <- rep(0,M)
  true_b2 <- rep(0,M)
  
  Jstat <- rep(0,M)
  
  for(i in 1:M){
    
    data <- dgp_clust(N = N,                  # Number of observations
                      G = G,                     # Number of clusters
                      p = 3,                      # Number of parameters in the parameter vector
                      endog = 2,                  # Number of endogenous regressors (first 'endog' columns of X)
                      serialcor = serialcor,              # (unused here)
                      clustersd = clustersd,              # SD of cluster common component in X
                      errors_var = errors_var,             # Var(idiosyncratic outcome error)
                      l = 4,                      # Number of instruments
                      rho_Z = rho_Z,                # Cov parameter for the idiosyncratic instrument component
                      sigmaZ = sigmaZ,                 # Variance of the idiosyncratic instrument component
                      alpha = 0.5,
                      rho_xe = 0.6,               # Corr between cluster shocks v_g (X) and u_g (Y)
                      indivsd = 1,                # SD of idiosyncratic X noise
                      sigmaZg = 1,                # Var of cluster-level instrument component Z_g
                      target_R2 = 0.4,            # Target first-stage R^2 for each endogenous regressor
                      sigma_u = 1.0,              # Var(cluster outcome shock u_g)
                      sigma_v = 1.0,              # Var(cluster regressor shock v_g)
                      seed = NULL,
                      beta = NULL                 # true beta; if NULL draw N(0,1)
                      )
    
    true_b1 <- data$beta[1]
    true_b2 <- data$beta[2]
    
    data <- data$data
    
    estimates <- gmm(Y ~ X1 + X2 + X3, 
                     ~ Z1 + Z2 + Z3 + Z4 + X3, 
                     data = data, 
                     type = "twoStep")
    
    coeffs <- estimates$coefficients
    
    vcov <- estimates[["vcov"]]
    
    lower_bound_int <- coeffs[[1]] - (z * sqrt(vcov[1,1])  / sqrt(N))
    upper_bound_int <- coeffs[[1]] + (z * sqrt(vcov[1,1])  / sqrt(N))
    coverage_int[i] <- (0.5>=lower_bound_int)*(0.5<=upper_bound_int)  
    
    lower_bound_b1 <- coeffs[[2]] - (z * sqrt(vcov[2,2])  / sqrt(N))
    upper_bound_b1 <- coeffs[[2]] + (z * sqrt(vcov[2,2])  / sqrt(N))
    coverage_b1[i] <- (true_b1>=lower_bound_b1)*(true_b1<=upper_bound_b1)  
    
    lower_bound_b2 <- coeffs[[3]] - (z * sqrt(vcov[3,3])  / sqrt(N))
    upper_bound_b2 <- coeffs[[3]] + (z * sqrt(vcov[3,3])  / sqrt(N))
    coverage_b2[i] <- (true_b2>=lower_bound_b2)*(true_b2<=upper_bound_b2)  
    
    Jstat[i] <- specTest(estimates)$test[1]
    
    bias_intercept[i] <- coeffs[[1]] - 0.5
    bias_b1[i] <- coeffs[[2]] - true_b1
    bias_b2[i] <- coeffs[[3]] - true_b2
    
  }
  
  return(list(biasint = mean(bias_intercept,na.rm = TRUE), 
              bias_b1 = mean(bias_b1, na.rm = TRUE),
              bias_b2 = mean(bias_b2, na.rm = TRUE),
              coverage_int = mean(coverage_int, na.rm = TRUE),
              coverage_b1 = mean(coverage_b1, na.rm = TRUE), 
              coverage_b2 = mean(coverage_b2, na.rm = TRUE), 
              J = Jstat
  ))
  
}

simulation(M=500, N = 500)

valid_inference <- function(M = 500, 
                            N = 1000,                  #Number of observations
                            G = 25,                     #Number of clusters
                            serialcor = 0.7,              #To allow for serial autocorrelation
                            clustersd = 1,              #Variance within each cluster (covariates)
                            errors_var = 1,             #Variance within each cluster (error terms)
                            rho_Z = 0.3,                #Parameter for the covariance of the instruments
                            sigmaZ = 1
                            ,level = 0.05){
  
  J_tilda <- rep(0,M)
  
  for(i in 1:M){
    
    data <- dgp_clust(N = N,                  # Number of observations
                      G = G,                     # Number of clusters
                      p = 3,                      # Number of parameters in the parameter vector
                      endog = 2,                  # Number of endogenous regressors (first 'endog' columns of X)
                      serialcor = serialcor,              # (unused here)
                      clustersd = clustersd,              # SD of cluster common component in X
                      errors_var = errors_var,             # Var(idiosyncratic outcome error)
                      l = 4,                      # Number of instruments
                      rho_Z = rho_Z,                # Cov parameter for the idiosyncratic instrument component
                      sigmaZ = sigmaZ,                 # Variance of the idiosyncratic instrument component
                      alpha = 0.5,
                      rho_xe = 0.6,               # Corr between cluster shocks v_g (X) and u_g (Y)
                      indivsd = 1,                # SD of idiosyncratic X noise
                      sigmaZg = 1,                # Var of cluster-level instrument component Z_g
                      target_R2 = 0.4,            # Target first-stage R^2 for each endogenous regressor
                      sigma_u = 1.0,              # Var(cluster outcome shock u_g)
                      sigma_v = 1.0,              # Var(cluster regressor shock v_g)
                      seed = NULL,
                      beta = NULL                 # true beta; if NULL draw N(0,1)
    )
    
    data <- data$data
    
    fs_gmm <- gmm(Y ~ X1 + X2 + X3, 
                  ~ Z1 + Z2 + Z3 + Z4 + X3, 
                  data = data, 
                  type = "twoStep")
    
    theta_hat <- fs_gmm$coefficients
    
    CCE <- CCE_Hwang(N = N, 
                     G = G, 
                     data = data, 
                     theta = theta_hat)
    
    J_init <- N * t(g(data = data, theta = theta_hat)) * solve(CCE) * g(data = data, theta = theta_hat)
    
    J_tilda[i] <- ((G - 1)/ 1) * J_init / (G - J_init)
  }
  
return(J_tilda)
  
}

valid_inference()
