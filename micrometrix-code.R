set.seed(1234)

options(error = NULL) #Avoids having the (annoying) R debugging tool at each error 
options(rlang_backtrace_on_error = "none")

install.packages(c("MASS", "parallel", "dplyr", "gmm", "Matrix"))
library(MASS)
library(parallel)
library(dplyr)
library(gmm)
library(Matrix)
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
  endog1 <- runif(1, 0, 0.4) 
  endog2 <- runif(1, 0, 0.4)
  
  ins_1 <- rnorm(l)
  
  X[,1] <- X[,1] +  Z %*% ins_1 + nu 
  
  ins_2 <- rnorm(l)
  
  X[,2] <- X[,2] + Z %*% ins_2 + eta
  
  # We generate the outcomes 
  
  ## Coefficients are randomly draw
  
  beta <- rnorm(p)
  
  Y <- X%*%beta + epsilon
  
  #We now generate our l instruments : 
  
  data <- data.frame(id = as.data.frame(id),Y = as.data.frame(Y), X = as.data.frame(X), Z = as.data.frame(Z))
  
  return(list(data = data,beta = beta))
  
}

testdata <- dgp_clust()
beta <- testdata$beta
testdata <- testdata$data

eps <- testdata$V1 - cbind(testdata$X.V1, testdata$X.V2, testdata$X.V3) %*% beta

mean(t(cbind(testdata$Z.V1, testdata$Z.V2, testdata$Z.V3, testdata$Z.V4)) %*% (testdata$V1 - cbind(testdata$X.V1, testdata$X.V2, testdata$X.V3) %*% beta))

cor(testdata$Z.V4, eps)




g <- function(beta, x, Z){
  
  g <- t(cbind(x$Z.V1, x$Z.V2, x$Z.V3, x$Z.V4)) %*% (x$V1 - cbind(x$X.V1, x$X.V2, x$X.V3) %*% beta)
  
}

estimates <- gmm(g, x = testdata, t0 = c(0,0), type = "twoStep")
