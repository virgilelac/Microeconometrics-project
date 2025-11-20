source("https://raw.githubusercontent.com/virgilelac/Microeconometrics-project/refs/heads/main/functions_sim.R")

M = 1000
n = 500
G = 25
L = 3
theta = 1
delta = 0.1
pi_vec = c(1.5, -2, 1.8)
sigma_a = 1
sigma_b = 1
sigma_eps = 1
sigma_eta = 1
level = 0.05

sim <- simulation(M,
                  n, G,
                  L, theta, delta, pi_vec,
                  sigma_a, sigma_b,sigma_eps, sigma_eta, level)

J_graph(J = sim$j)

ks.test(sim$j, "pchisq", df = 2)

sum(sim$j < qchisq(0.95, df = 2))/1000

