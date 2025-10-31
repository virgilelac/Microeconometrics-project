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
  u <- drop(data$Y - X %*% theta)                       
  
  (Z*u)                                                   
}

CCE_Hwang <- function(N, G, data, theta){
  
  ## We compute g (expressed as f in Hwang 2021) for each individual 
  
  normalize <- 1/sqrt((1/G) * sum(as.integer(table(data$id))))
  
  g_hat <- as.data.frame(g(data, theta))
  g_hat$id <- data$id
  sum_i <- g_hat %>%
    group_by(id) %>%
    summarize(across(everything(), sum)) %>%
    mutate(across(-id, ~ . * normalize)) %>% 
    select(-id)
  
  sum_i <- as.matrix(sum_i)
  
  sum_i_t <- lapply(1:nrow(sum_i), function(i) {
    sum_vector <- as.vector(sum_i[1,])
    tcrossprod(sum_vector)
  })
  
  CCE <- Reduce("+", sum_i_t)/G
  
  return(as.matrix(CCE))
  
}

#################################################################################################
############################ HIGHER LEVEL FUNCIONS ##############################################
#################################################################################################

dgp_clust <-function( 
    N = 10000,                  #Number of observations
    G = 25,                     #Number of clusters
    p = 3,                      #Number of parameters in the parameter vector
    endog = 3,                  #Number of endogenous regressors
    serialcor = 0,              #To allow for serial autocorrelation
    clustersd = 1,              #Variance within each cluster (covariates)
    errors_var = 1,             #Variance within each cluster (error terms)
    l = 4,                      #Number of instruments
    rho_Z = 0.3,                #Parameter for the covariance of the instruments
    sigmaZ = 1                  #Variance of the instruments
){
  
  if(l<=p){stop("Number of instruments should exceed the number of parameters")}
  
  ### Instruments
  
  S_ins <- varcov_ins(l, rho_Z, sigmaZ)
  
  Z <- mvrnorm(n = N, mu = rep(0, l), Sigma = S_ins)
  
  ### Group structure
  
  # To allow for unbalanced clusters, keeping approximately same sizes (as in Hwang 2021)
  
  id <- groupmaker(N,G)
  
  #We generate covariates and error terms 
  
  S_endog <- varcov_ins(endog, 0.6, 1)
  
  mu_endog <- rep(0, endog)
  
  errors <- mvrnorm(N, mu = mu_endog, Sigma = S_endog)
  
  #We generate arbitrary correlation bewteen the endogenous regressors and the noise term 
  
  epsilon <- as.vector(errors[,1])
  nu <- as.vector(errors[,2])
  eta <- as.vector(errors[,3])
  
  X <- covariates_gen(id, N, p, G, clustersd)
  
  #We randomly assign the degree to which the endogenous regressors are endognenous
  #The endogenous regressors are arbitrarily set to be the 2 first regressors
  
  ins_1 <- rnorm(l)
  
  X[,1] <- X[,1] +  Z %*% ins_1 + nu 
  
  ins_2 <- rbeta(l, 0.7, 0.5)
  
  X[,2] <- X[,2] + Z %*% ins_2 + eta
  
  # We generate the outcomes 
  
  ## Coefficients are randomly draw
  
  beta <- rnorm(p)
  
  Y <- X%*%beta + epsilon
  
  #We now generate our l instruments : 
  
  colnames(X) <- paste0("X", seq_len(ncol(X)))
  colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
  data <- data.frame(id = id, Y, X, Z)
  
  return(list(data = data,beta = beta))
  
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
    
    data <- dgp_clust(N = N,                  #Number of observations
                      G = G,                     #Number of clusters
                      p = 3,                      #Number of parameters in the parameter vector
                      endog = 3,                  #Number of endogenous regressors
                      serialcor = serialcor,              #To allow for serial autocorrelation
                      clustersd = clustersd,              #Variance within each cluster (covariates)
                      errors_var = errors_var,             #Variance within each cluster (error terms)
                      l = 4,                      #Number of instruments
                      rho_Z = rho_Z,                #Parameter for the covariance of the instruments
                      sigmaZ = sigmaZ                 # true beta; if NULL draw N(0,1)
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
    
    data <- dgp_clust(N = N,                  #Number of observations
                      G = G,                     #Number of clusters
                      p = 3,                      #Number of parameters in the parameter vector
                      endog = 3,                  #Number of endogenous regressors
                      serialcor = serialcor,              #To allow for serial autocorrelation
                      clustersd = clustersd,              #Variance within each cluster (covariates)
                      errors_var = errors_var,             #Variance within each cluster (error terms)
                      l = 4,                      #Number of instruments
                      rho_Z = rho_Z,                #Parameter for the covariance of the instruments
                      sigmaZ = sigmaZ                 # true beta; if NULL draw N(0,1)
    )
    
    data <- data$data
    
    fs_gmm <- gmm(Y ~ X1 + X2 + X3, 
                  ~0 + Z1 + Z2 + Z3 + Z4 + X3, 
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

dataset <- dgp_clust(N = 10000)
theta <- dataset$beta
data <- dataset$data

fs_gmm <- gmm(Y ~ X1 + X2 + X3, 
              ~0 + Z1 + Z2 + Z3 + Z4 + X3, 
              data = data, 
              type = "twoStep")

theta_hat <- fs_gmm$coefficients
theta <- c(0.5,dataset$beta)

CCE <- CCE_Hwang(N = 100000, 
                 G = 25, 
                 data = data, 
                 theta = theta_hat)

g(dataset, theta_hat)


sum(split_g_hat[["1"]])

