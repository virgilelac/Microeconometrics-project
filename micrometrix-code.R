options(error = NULL) #Avoids having the (annoying) R debugging tool at each error 
options(rlang_backtrace_on_error = "none")

install.packages(c("MASS", "parallel", "dplyr", "gmm"))

library(MASS)
library(parallel)
library(dplyr)
library(gmm)
#### Replication code for the microeconometrics project. 

### DGP function to create clustered data 

dgp_clust <- function(N = 1000, theta, intercept, clustering = FALSE){ 
  
  ## Creation of the instruments 
  
  ins <- mvrnorm(n = N, mu = c(0, 0), Sigma = matrix(c(2, 0.8, 0.8, 3), ncol = 2))
  
  Z1 <- as.vector(ins[,1])
  Z2 <- as.vector(ins[,2])
  
  ## creation of the endogenous variable
  
  eps <- rnorm(N, mean = 0, sd = 2)
  
  rho <- runif(1, 0,1)
  
  u <- rnorm(N, mean = 0, sd = 2)
  
  u <- rho * eps + sqrt((1 - rho^2)) * u
  
  X <- 1.3 * Z1 + 0.8 * Z2 + u
  
  if(clustering == TRUE){ 
    
    ## creation of the endogenous variable
    
    eps1 <- rnorm(N/2, mean = 0, sd = 1)
    
    rho <- runif(1, 0,1)
    
    u1 <- rnorm(N/2, mean = 0, sd = 1)
    
    u1 <- rho * eps1 + sqrt((1 - rho^2)) * u1
    
    X1 <- 1.3 * Z1[1:(N/2)] + 0.8 * Z2[1:(N/2)] + u1
    
    eps2 <- rnorm(N/2, mean = 0, sd = 4)
    
    u2 <- rnorm(N/2,mean = 0, sd = 1)
    
    u2 <- rho*eps2 + sqrt((1 - rho^2)) + u2 
    
    X2 <- 1.3 * Z1[((N/2) + 1):N] + 0.8 * Z2[((N/2) + 1):N] + u2
    
    X <- as.vector(rbind(X1, X2))
    
    eps <- as.vector(rbind(eps1, eps2))
    
    }
  
  ## creation of the outcomes 
  
  Y <- intercept + theta * X + eps
  
  data <- data.frame(Y = Y, X = X, Z1 = Z1, Z2 = Z2)
  
  return(list(data = data, theta = theta, intercept = intercept))
  
}

g <- function(theta, x){
  
  g <- cbind(1, x$Z1, x$Z2) * (x$Y - theta[1] - theta[2] * x$X)
  
}

simulation <- function(M, N, clustering = FALSE, level = 0.05){
  
  theta_sim <- rnorm(1, mean = 2, sd = 1)
  
  intercept_sim <- rnorm(1,mean = 0.5, sd = 0.25)
  
  theta <- rep(0,M)
  intercept <- rep(0,M)
  se_inter <- rep(0,M)
  se_t <- rep(0,M)
  coverage <- rep(0,M)
  
  z <- qnorm(1-level/2)
  
  Jstat <- rep(0,M)
  
  for(i in 1:M){
    
    data <- dgp_clust(N = N, theta = theta_sim, intercept = intercept_sim, clustering = clustering)$data
    
    estimates <- gmm(g, x = data, t0 = c(0,0), type = "twoStep")
    
    coeffs <- estimates$coefficients
    
    intercept[i] <- coeffs[[1]]
    theta[i] <- coeffs[[2]]
  
    vcov <- estimates[["vcov"]]
    
    se_t[i] <- sqrt(vcov[2,2]) 
    
    se_inter[i] <- sqrt(vcov[1,1])
    
    lower_bound <- coeffs[[2]] - (z * sqrt(vcov[2,2])  / sqrt(N))
    upper_bound <- coeffs[[2]] + (z * sqrt(vcov[2,2])  / sqrt(N))
    
    coverage[i] <- (theta_sim>=lower_bound)*(theta_sim<=upper_bound)  
    
    Jstat[i] <- specTest(estimates)$test[1]
    
  }
  
  return(list(intercept_true = intercept_sim, theta_true = theta_sim, intercept = mean(intercept), theta = mean(theta), coverage = mean(coverage), Jstats = Jstat))
  
}

graph <- function(statistics, distribution = "CDF"){ 
  
  J  <- statistics
  df <- 1 
  
  if(distribution == "CDF"){
    
    plot <- plot(ecdf(J))
    curve(pchisq(x, df = df), add = TRUE, lwd = 2)
    legend("bottomright", c("Empirical CDF", "Theoretical CDF"),
           lwd = c(1,2), col = c("black","black"), bty = "n")
    
  }
  
  if(distribution == "PDF"){
    
    d <- density(J)
    plot(d)          # empirical curve
    curve(dchisq(x,df = df), add = TRUE, lwd = 2)  # theoretical curve
    legend("topright", c("Empirical KDE", "Theoretical PDF"),
           lwd = c(1,2), bty = "n")
    
  }
  
  }


nonclustered <- simulation(M = 1000, N = 1000, clustering = FALSE, level = 0.05)

clustered <- simulation(M = 1000, N = 1000, clustering = TRUE, level = 0.05)

graph(clustered$Jstats)

graph(nonclustered$Jstats)
