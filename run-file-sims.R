source("https://raw.githubusercontent.com/virgilelac/Microeconometrics-project/refs/heads/main/functions_sim.R")

#### Wrongful inference: consider no clustering when there is clustering. 

table <- sim_table(M = 500,n_vec = c(3000, 6000, 9000), G_vec = c(20, 30), 
          L = 3,
          theta = 1,
          delta = 0.3,
          pi_vec = c(1.5, -2, 1.8),
          sigma_a = 1,
          sigma_b = 1,
          sigma_eps = 1,
          sigma_eta = 1,
          rho = 0.1,
          level = 0.05)

latex_table(table,
            label   = "tab:sec6",
            caption = "Simulation results: Hwang's case versus neglecting cluster structure")

### Random little checks


lala <- ks.test(sim$j, "pchisq", df = 2)$p.value

J_graph(J = sim$j_Hwang)

ks.test(sim$j_Hwang, "pf", df1 = 2, df2 = 3
        )


