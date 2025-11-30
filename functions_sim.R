set.seed(1234)

options(error = NULL) #Avoids having the (annoying) R debugging tool at each error 
options(rlang_backtrace_on_error = "none")

packs <- c("MASS", "parallel", "dplyr", "gmm", "Matrix", "future", "future.apply", "torch")
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
library(torch)

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

CCE_Hwang_centered <- function(N, G, data, theta){

  g_hat <- as.data.frame(g(data, theta))
  g_hat$id <- data$id

  normalize <- 1 / sqrt((1 / G) * sum(as.integer(table(data$id))))

  sum_i <- g_hat %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(
      dplyr::across(
        .cols = dplyr::everything(), 
        .fns  = ~ sum(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )
  
  H <- as.matrix(sum_i[ , -1, drop = FALSE])

  H_centered <- normalize * scale(H, center = TRUE, scale = FALSE) 

  CCE_centered <- crossprod(H_centered) / (G)
  
  return(as.matrix(CCE_centered))
}

### Main functions

dgp_clust <- function(n = 10000, G = 5,
                      L = 3, theta = 1, delta = 0.1, pi_vec = c(0.8, 0.5, 0.3),
                      sigma_a = 1, sigma_b = 1,sigma_eps = 1, sigma_eta = 1, rho = 0){
  
  ###Regularity checks
  
  if(length(pi_vec) != L){stop("pi must be of length L")}
  if(n%%G != 0){stop("n%%G must be close to 0")}
  
  ### Observations per cluster 
  
  m <- n/G
  
  ### Assignation to clusters 
  
  id <- rep(1:G, each = m)
  
  # ### Cluster-level shocks

  a_g <- rnorm(G, 0, sigma_a)     # a_g ~ N(0, sigma_a^2)
  b_g <- rnorm(G, 0, sigma_b)     # b_g ~ N(0, sigma_b^2)
  
  # ### Expand cluster shocks to observation level
  
  a <- a_g[id]
  b <- b_g[id]

  ### Draw individual shocks
  eta <- rnorm(n, 0, sigma_eta)

  eps <- rnorm(n, 0, sigma_eps)
  
  ### Draw instruments Z (n x L)
  zeta_g   <- matrix(rnorm(G * L), nrow = G, ncol = L)
  zeta_eps <- matrix(rnorm(n * L), nrow = n, ncol = L)
  Z <- rho * zeta_g[id, ] + sqrt(1 - rho^2) * zeta_eps 
  
  ### First-stage error: v = b + eta
  #v <- eta
  ### Structural error: u = delta*v + a + eps
  v <- b +  eta
  u <- delta * v + a +  eps
  
  ### Generate X and Y
  x <- Z %*% pi_vec + v 
  y <- x * theta + u                   # Y = theta*X + u
  
  #Create data frame 
  data <- data.frame(
    id = id,
    y = as.numeric(y),
    x = as.numeric(x),
    Z1 = Z[,1],
    Z2 = Z[,2],
    Z3 = Z[,3]
  )
  
  return(data)
  
}

simulation <- function(M = 1000, 
                       n = 50000, G = 25,
                       L = 3, theta = 1, delta = 0.1, pi_vec = c(0.8, 0.5, 0.3),
                       sigma_a = 1, sigma_b = 1,sigma_eps = 1, sigma_eta = 1, rho = 0,level = 0.05){
  
  theta_vec <- rep(0,M)
  theta_Hwang <- rep(0,M)
  bias <- rep(0,M)
  bias_Hwang <- rep(0,M)
  cov_j <- rep(0,M)
  j_vec <- rep(0,M)
  j_Hwang <- rep(0,M)
  cov_Hwang <- rep(0,M)
  
  true_theta <- theta
  q <- L - length(theta)
  
  m <- qchisq(1 - level, df = q)
  f <- qf(1 - level, df1 = q, df2 = G - q)
  
  for (i in 1:M) {
    
    
    data <- dgp_clust(n, G,
                      L, theta, delta, pi_vec,
                      sigma_a, sigma_b,sigma_eps, sigma_eta, rho)
    
    model <- gmm(y ~ 0 + x, 
                 ~ 0 + Z1 + Z2 + Z3, 
                 data = data, 
                 weightsMatrix = diag(L))
    
    eff <- gmm(y ~ 0 + x, 
               ~ 0 + Z1 + Z2 + Z3, 
               data = data, 
               method = "twoStep")
    
    estimates <- model$coefficients
    
    vcov <- solve((CCE_Hwang_centered(n,G,data,estimates)))
    
    twostep <- gmm(y ~ 0 + x,
                   ~ 0 + Z1 + Z2 + Z3, 
                   data = data, 
                   weightsMatrix = vcov
    )
    
    estimates <- eff$coefficients
    estimates_Hwang <- twostep$coefficients
 
    j <- specTest(eff)$test[1]
    j_vec[i] <- j
    
    g_mat  <- g(data, coef(twostep))              
    g_bar  <- colMeans(g_mat)             
    J    <- as.numeric(n * t(g_bar) %*% vcov %*% g_bar)
    
    j_tilda <- ((G - q)/(G*q)) * J
    j_Hwang[i] <- j_tilda
    
    theta_Hwang[i] <- estimates_Hwang[1]
    theta_vec[i] <- estimates[1]

    bias[i] <- estimates[1] - true_theta[1]
    bias_Hwang[i] <- estimates_Hwang[1] - true_theta[1]
  
    cov_j[i] <- specTest(eff)$test[1] < m
    cov_Hwang[i] <- j_tilda < f
    
  }
  
  RMSE      <- sqrt(mean((theta_vec   - true_theta)^2))
  RMSE_hwang <- sqrt(mean((theta_Hwang - true_theta)^2))
  
  results <- list(
    theta = mean(theta_vec, na.rm = TRUE),
    theta_hwang = mean(theta_Hwang, na.rm = TRUE),
    bias = mean(bias, na.rm = TRUE),
    bias_Hwang = mean(bias_Hwang, na.rm = TRUE),
    RMSE = RMSE,
    RMSE_hwang = RMSE_hwang,
    j        = j_vec, 
    j_Hwang = j_Hwang,
    cov_j = mean(cov_j, na.rm = TRUE),
    cov_J_Hw = mean(cov_Hwang, na.rm = TRUE)
  )
  
  return(results)
}

sim_table <- function(M = 1000,n_vec = c(100, 500, 1000), G_vec = c(4, 5), 
                      L = 3, theta = 1, delta = 0.1, pi_vec = c(0.8, 0.5, 0.3),
                      sigma_a = 1, sigma_b = 1,sigma_eps = 1, sigma_eta = 1, rho = 0.1,level = 0.05){
  
  g_vals <- as.vector(G_vec)
  n_vals <- as.vector(n_vec)
  
  param_grid <- expand.grid(
    g = g_vals,
    n = n_vals
  )
  
  cl <- makeCluster(detectCores()- 2, type = "FORK")
  
    res <- pblapply(1:nrow(param_grid), cl = cl, function(i){
                      
      row <- as.numeric(param_grid[i, ])
      
      n <- row[2]
      G <- row[1]
      sim <- simulation(M,
                       n, G,
                      L, theta, delta, pi_vec,
                      sigma_a, sigma_b,sigma_eps, sigma_eta, rho,level)
      
      ks_std <-  ks.test(sim$j, "pchisq", df = L - length(theta))$p.value > 0.05
      ks_Hw <- ks.test(sim$j_Hwang, "pf", df1 = 2, df2 = G - (L - length(theta)))$p.value > 0.05
                      
      res <- data.frame(         n = row[2],
                                 g = row[1], 
                                 theta = sim$theta_hwang,
                                 bias = sim$bias_Hwang, 
                                 RMSE = sim$RMSE, 
                                 Std_cov = sim$cov_j, 
                                 Hwang_cov = sim$cov_J_Hw, 
                                 ks_std = ks_std,
                                 ks_Hw = ks_Hw
                                  )
  })
      
    stopCluster(cl)
    
    results_df <- do.call(rbind, res)
    
    return(results_df)
    
}


J_graph <- function(J){
  
  df <- data.frame(J = J)
  
  ggplot(df, aes(x = J)) +
    geom_density() +
    labs(x = "J", y = "Density") +
    #coord_cartesian(xlim = c(-20, 100)) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
}

latex_table <- function(results, 
                        label   = "tab:sim_results", 
                        caption = "Simulation results") {
  # results: data.frame from sim_table()
  # Columns expected: n, g, bias, RMSE, Std_cov, Hwang_cov, ks_std, ks_Hw
  
  library(dplyr)
  
  # Ensure proper ordering
  results_df <- results %>%
    arrange(g, n)
  
  g_vals <- sort(unique(results_df$g))
  n_vals <- sort(unique(results_df$n))
  
  cat("\\begin{table}[H]\n")
  cat("\\centering\n")
  cat(paste0("\\caption{", caption, "}\n"))
  cat(paste0("\\label{", label, "}\n"))
  cat("\\resizebox{0.7\\textwidth}{!}{%\n")
  
  # Header
  cat("\\begin{tabular}{ll", paste(rep("r", length(n_vals)), collapse = ""), "}\n", sep = "")
  cat("\\toprule\n")
  
  cat("$G$ & Statistic")
  for (n_val in n_vals) {
    cat(paste0(" & $n=", n_val, "$"))
  }
  cat(" \\\\\n")
  cat("\\midrule\n")
  
  # Helper for ks-test labels
  ks_label <- function(x) {
    if (is.na(x)) {
      return("---")
    }
    if (isTRUE(x)) {
      return("No rejection")
    } else {
      return("Rejection")
    }
  }
  
  for (g_val in g_vals) {
    # Order of statistics in each G-block
    stat_order <- c("Bias", 
                    "RMSE", 
                    "Uncorrected cov.", 
                    "Uncor. ks test", 
                    "Corrected cov.", 
                    "Corr. ks test")
    
    for (k in seq_along(stat_order)) {
      stat_name <- stat_order[k]
      
      # First row of block: show G, others use "-"
      if (k == 1) {
        cat(g_val, " & ", stat_name, sep = "")
      } else {
        cat("- & ", stat_name, sep = "")
      }
      
      for (n_val in n_vals) {
        val <- results_df %>% 
          filter(g == g_val, n == n_val)
        
        if (nrow(val) == 0) {
          cat(" & ---")
        } else {
          if (stat_name == "Bias") {
            cat(" & ", sprintf("%.3f", val$bias), sep = "")
          } else if (stat_name == "RMSE") {
            cat(" & ", sprintf("%.3f", val$RMSE), sep = "")
          } else if (stat_name == "Uncorrected cov.") {
            cat(" & ", sprintf("%.3f", val$Std_cov), sep = "")
          } else if (stat_name == "Corrected cov.") {
            cat(" & ", sprintf("%.3f", val$Hwang_cov), sep = "")
          } else if (stat_name == "Uncor. ks test") {
            cat(" & ", ks_label(val$ks_std), sep = "")
          } else if (stat_name == "Corr. ks test") {
            cat(" & ", ks_label(val$ks_Hw), sep = "")
          }
        }
      }
      cat(" \\\\\n")
    }
    
    cat("\\midrule\n")
  }
  
  cat("\\bottomrule\n")
  cat("\\end{tabular}%\n")
  cat("}\n")
  cat("\\end{table}\n")
}
