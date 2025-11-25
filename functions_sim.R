set.seed(1234)

options(error = NULL) #Avoids having the (annoying) R debugging tool at each error 
options(rlang_backtrace_on_error = "none")

packs <- c("MASS", "parallel", "dplyr", "gmm", "Matrix", "future", "future.apply")
inst <- packs[!packs %in% rownames(installed.packages())]
if (length(inst) > 0) install.packages(inst)

library(MASS)
library(parallel)
library(dplyr)
library(gmm)
library(Matrix)
library(ggplot2)
library(future)
library(future.apply)
library(pbapply)

### Hwang's CCE 

g <- function(data,theta) {
  X <- model.matrix(~ 0 + x, data = data)              
  Z <- model.matrix(~ 0 + Z1 + Z2 + Z3, data = data)    
  u <- drop(data$y - X %*% theta)                       
  
  (Z*u)                                                   
}

CCE_Hwang <- function(N, G, data, theta){
  
  ## We compute g (expressed as f in Hwang 2021) for each individual 
  
  normalize <- 1/sqrt((1/G) * sum(as.integer(table(data$id))))
  
  g_hat <- as.data.frame(g(data, theta))
  g_hat$id <- data$id
  sum_i <- g_hat %>%
    group_by(id) %>%
    summarize(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
    mutate(across(-id, ~ . * normalize)) %>% 
    select(-id)
  
  sum_i <- as.matrix(sum_i)
  
  sum_i_t <- lapply(1:nrow(sum_i), function(i) {
    sum_vector <- as.vector(sum_i[i,])
    tcrossprod(sum_vector)
  })
  
  CCE <- Reduce("+", sum_i_t)/G
  
  return(as.matrix(CCE))
  
}

### Main functions

dgp_clust <- function(n = 50000, G = 25,
                      L = 3, theta = 1, delta = 0.1, pi_vec = c(0.8, 0.5, 0.3),
                      sigma_a = 1, sigma_b = 1,sigma_eps = 1, sigma_eta = 1){
  
  ###Regularity checks
  
  if(length(pi_vec) != L){stop("pi must be of length L")}
  if(n%%G != 0){stop("n%%G must be close to 0")}
  
  ### Observations per cluster 
  
  m <- n/G
  
  ### Assignation to clusters 
  
  id <- rep(1:G, each = m)
  
  ### Cluster-level shocks
  
  a_g <- rnorm(G, 0, sigma_a)     # a_g ~ N(0, sigma_a^2)
  b_g <- rnorm(G, 0, sigma_b)     # b_g ~ N(0, sigma_b^2)
  
  ### Expand cluster shocks to observation level
  a <- a_g[id]
  b <- b_g[id]
  
  ### Draw individual shocks
  eta <- rnorm(n, 0, sigma_eta)
  eps <- rnorm(n, 0, sigma_eps)
  
  ### Draw instruments Z (n x L)
  Z <- matrix(rnorm(n * L), n, L)  #NOTE Z a standard normal - could change this later? 
  
  ### First-stage error: v = b + eta
  v <- b + eta
  #v <- eta
  ### Structural error: u = delta*v + a + eps
  u <- delta * v + a + eps
  
  ### Generate X and Y
  x <- Z %*% pi_vec + v                # X = Z*pi + v
  y <- x * theta + u                   # Y = theta*X + u
  
  #Create data frame 
  dat <- data.frame(
    id = id,
    y = as.numeric(y),
    x = as.numeric(x),
    Z1 = Z[,1],
    Z2 = Z[,2],
    Z3 = Z[,3]
  )
}

simulation <- function(M = 10000, 
                       n = 50000, G = 25,
                       L = 3, theta = 1, delta = 0.1, pi_vec = c(0.8, 0.5, 0.3),
                       sigma_a = 1, sigma_b = 1,sigma_eps = 1, sigma_eta = 1, level = 0.05){
  
  step <- function(n, G,
                   L, theta, delta, pi_vec,
                   sigma_a, sigma_b,sigma_eps, sigma_eta, level){
    
    z <- qnorm(1-level/2) 
    
    data <- dgp_clust(n, G,
                      L, theta, delta, pi_vec,
                      sigma_a, sigma_b,sigma_eps, sigma_eta)
    
    model <- gmm(y ~0 + x, 
                 ~ 0 + Z1 + Z2 + Z3, 
                 data = data, 
                 type = "twoStep")
    
    true_theta <- theta
    
    estimates <- model$coefficients
    
    vcov <- model$vcov
    
    lower <- estimates[1] - z *sqrt(as.numeric(vcov[1]))
    upper <- estimates[1] + z *sqrt(as.numeric(vcov[1]))
    
    theta_est <- estimates[1]
    coverage <- (true_theta >= lower)*(true_theta<= upper)
    jstat <- specTest(model)$test[1]
    
    cvg <- jstat < 3.84
    
    return(list(theta = theta_est,
                coverage = coverage,
                j = jstat,
                cvg = cvg))
    
  }
  
  n = n
  G = G
  L = L
  theta = theta
  delta = delta
  pi_vec = pi_vec
  sigma_a = sigma_a
  sigma_b = sigma_b
  sigma_eps = sigma_eps
  sigma_eta = sigma_eta
  level = level
  
  cl <- makeCluster(min(20, detectCores() - 2), type = if(.Platform$OS.type == "windows"){"PSOCK"}else{"FORK"})
  
  ## Configuration for windows 
  
  if(.Platform$OS.type == "windows"){
    clusterEvalQ(cl, {library(MASS)
      library(parallel)
      library(dplyr)
      library(gmm)
      library(Matrix)
    })
    clusterExport(cl, varlist = c("dgp_clust", "step"),envir = environment())} 
  
  sim <- pblapply(1:M, cl = cl, function(i) {set.seed(1234+i) ## Given the type of cluster we use for parallelization, we need this 
    step(n, 
         G,
         L,
         theta, 
         delta,
         pi_vec,
         sigma_a, 
         sigma_b,
         sigma_eps, 
         sigma_eta, 
         level)
  })
  
  
  theta_vec <- sapply(sim, `[[`, "theta")
  coverage_vec <- sapply(sim, `[[`, "coverage")
  j_vec <- sapply(sim, `[[`, "j")
  cov_j <- sapply(sim, `[[`, "cvg")
  
  results <- list(
    theta    = mean(theta_vec, na.rm = TRUE),
    coverage = mean(coverage_vec, na.rm = TRUE),
    j        = j_vec, 
    cov_j = mean(cov_j, na.rm = TRUE)
  )
  
  return(results)
}

J_graph <- function(J){
  
  df <- data.frame(J = J)
  
  ggplot(df, aes(x = J)) +
    geom_density() +
    labs(x = "J", y = "Density") +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
}

good_inference <- function(M = 10000, 
                           n = 50000, G = 25,
                           L = 3, theta = 1, delta = 0.1, pi_vec = c(0.8, 0.5, 0.3),
                           sigma_a = 1, sigma_b = 1,sigma_eps = 1, sigma_eta = 1, level = 0.05){
  
  step <- function(n, G,
                   L, theta, delta, pi_vec,
                   sigma_a, sigma_b,sigma_eps, sigma_eta, level){
    
    q <- L - length(theta)
    data <- dgp_clust(n, G,
                      L, theta, delta, pi_vec,
                      sigma_a, sigma_b,sigma_eps, sigma_eta)
    
    model <- gmm(y ~ 0 + x, 
                 ~ 0 + Z1 + Z2 + Z3, 
                 data = data, 
                 type = "twoStep")
    
    true_theta <- theta
    
    estimates <- model$coefficients
    
    vcov <- CCE_Hwang(N,G,data,estimates)
    
    twostep <- gmm(y ~ 0 + x,
                   ~ 0 + Z1 + Z2 + Z3, 
                   data = data, 
                   weightsMatrix = vcov
    )
    
    estimates <- twostep$coefficients
    
    theta_est <- estimates[1]
    
    jstat <- specTest(model)$test[1]
    
    j_tilda <- ((G - q)/q) * (q*jstat/(G - q*jstat))
    
    cvg <- j_tilda <= qf(1 - level, df1 = 2, df2 = G - 2)
    
    return(list(theta = theta_est,
                j = j_tilda,
                cvg = cvg))
    
  }
  
  n = n
  G = G
  L = L
  theta = theta
  delta = delta
  pi_vec = pi_vec
  sigma_a = sigma_a
  sigma_b = sigma_b
  sigma_eps = sigma_eps
  sigma_eta = sigma_eta
  level = level
  
  cl <- makeCluster(min(20, detectCores() - 2), type = if(.Platform$OS.type == "windows"){"PSOCK"}else{"FORK"})
  
  ## Configuration for windows 
  
  if(.Platform$OS.type == "windows"){
    clusterEvalQ(cl, {library(MASS)
      library(parallel)
      library(dplyr)
      library(gmm)
      library(Matrix)
    })
    clusterExport(cl, varlist = c("dgp_clust"),envir = environment())} 
  
  sim <- pblapply(1:M, cl = cl, function(i) {set.seed(1234+i) ## Given the type of cluster we use for parallelization, we need this 
    step(n, 
         G,
         L, 
         theta, 
         delta, 
         pi_vec,
         sigma_a, sigma_b,sigma_eps, sigma_eta, 
         level)})
  
  
  theta_vec <- sapply(sim, `[[`, "theta")
  j_vec <- sapply(sim, `[[`, "j")
  cov_j <- sapply(sim, `[[`, "cvg")
  
  results <- list(
    theta    = mean(theta_vec, na.rm = TRUE),
    j        = j_vec, 
    cov_j = mean(cov_j, na.rm = TRUE)
  )
  
  return(results)
  
}
