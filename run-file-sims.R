source("https://raw.githubusercontent.com/virgilelac/Microeconometrics-project/refs/heads/main/functions_sim.R")

#### Wrongful inference: consider no clustering when there is clustering. 
set.seed(1)
table_case1 <- sim_table(M = 500,n_vec = c(3000, 6000, 9000), G_vec = c(20, 30), 
          L = 3,
          theta = 1,
          delta = 0.3,
          pi_vec = c(1.5, -2, 1.8),
          sigma_a = 1,
          sigma_b = 1,
          sigma_eps = 1,
          sigma_eta = 1,
          rho = 0.1,
          level = 0.05, 
          method = "classic")

latex_table(table_case1,
            label   = "tab:sec6.1",
            caption = "Simulation results: Hwang's case versus neglecting cluster structure")

## Case 2 

table_case2 <- sim_table(M = 500,n_vec = c(3000, 6000, 9000), G_vec = c(20, 30), 
                         L = 3,
                         theta = 1,
                         delta = 0.3,
                         pi_vec = c(1.5, -2, 1.8),
                         sigma_a = 1,
                         sigma_b = 1,
                         sigma_eps = 1,
                         sigma_eta = 1,
                         rho = 0.1,
                         level = 0.05, 
                         method = "White")

latex_table(table_case2,
            label   = "tab:sec6.2",
            caption = "Simulation results: Hwang's case versus standard large-G asymptotics")


## The weak instruments problem 

table_weak <- sim_table(M = 500,n_vec = c(3000, 6000, 9000), G_vec = c(20, 30), 
                         L = 3,
                         theta = 1,
                         delta = 0.3,
                         pi_vec = c(0.005, -0.03, 0.001),
                         sigma_a = 1,
                         sigma_b = 1,
                         sigma_eps = 1,
                         sigma_eta = 1,
                         rho = 0.1,
                         level = 0.05, 
                         method = "White")

table_weak_difflev <- sim_table(M = 500,n_vec = c(3000, 6000, 9000), G_vec = c(20, 30), 
                        L = 3,
                        theta = 1,
                        delta = 0.3,
                        pi_vec = c(0.005, -0.03, 0.001),
                        sigma_a = 1,
                        sigma_b = 1,
                        sigma_eps = 1,
                        sigma_eta = 1,
                        rho = 0.1,
                        level = 0.1, 
                        method = "White")

latex_table(table_weak,
            label   = "tab:sec6.3",
            caption = "Simulation results: weak instruments")


### The heterogenous clusters problem

table_het <- sim_table(M = 500,n_vec = c(3000, 6000, 9000), G_vec = c(20, 30), 
                        L = 3,
                        theta = 1,
                        delta = 0.3,
                        pi_vec = c(1.5, -2, 1.8),
                        sigma_a = 1,
                        sigma_b = 1,
                        sigma_eps = 1,
                        sigma_eta = 1,
                        rho = 0.1,
                        level = 0.05, 
                        method = "White", 
                       heterogenous = TRUE)

table_het_difflev <- sim_table(M = 500,n_vec = c(3000, 6000, 9000), G_vec = c(10, 15), 
                       L = 3,
                       theta = 1,
                       delta = 0.3,
                       pi_vec = c(1.5, -2, 1.8),
                       sigma_a = 1,
                       sigma_b = 1,
                       sigma_eps = 1,
                       sigma_eta = 1,
                       rho = 0.1,
                       level = 0.05, 
                       method = "White", 
                       heterogenous = TRUE)

latex_table(table_het,
            label   = "tab:sec6.4",
            caption = "Simulation results: heterogenous clusters")

testhet <- simulation(M = 500,n = 10000, G = 25, 
           L = 3,
           theta = 1,
           delta = 0.3,
           pi_vec = c(1.5, -2, 1.8),
           sigma_a = 1,
           sigma_b = 1,
           sigma_eps = 1,
           sigma_eta = 1,
           rho = 0.1,
           level = 0.05, 
           method = "White", 
           heterogenous = TRUE)
