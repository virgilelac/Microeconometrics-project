source("https://raw.githubusercontent.com/virgilelac/Microeconometrics-project/refs/heads/main/functions_sim.R")

#### Wrongful inference: consider no clustering when there is clustering. 

M = 10000
n = 50000
G = 10
L = 3
theta = 1
delta = 0.1
pi_vec = c(1.5, -2, 1.8)
sigma_a = 1
sigma_b = 1
sigma_eps = 1
sigma_eta = 1
level = 0.05

sim_wrong <- simulation(M,
                  n, G,
                  L, theta, delta, pi_vec,
                  sigma_a, sigma_b,sigma_eps, sigma_eta, level)

J_graph(J = sim_wrong$j)

ks.test(sim$j, "pchisq", df = 2)


#### Good inference: assuming clustering (and fixed-G asymptotics) when it's the case

M = 10000
n = 50000
G = 10
L = 3
theta = 1
delta = 0.1
pi_vec = c(1.5, -2, 1.8)
sigma_a = 1
sigma_b = 1
sigma_eps = 1
sigma_eta = 1
level = 0.05

sim_good <- good_inference(M,
                       n, G,
                       L, theta, delta, pi_vec,
                       sigma_a, sigma_b,sigma_eps, sigma_eta, level)

J_graph(J = sim_good$j)

ks.test(sim_good$j, "pf", df1 = 2, df2 = G - 2)
