set.seed(1234)

M = 10000
n = 50000
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

simulation_bad_j_stat <- simulation(M,
                                    n, G,
                                    L, theta, delta, pi_vec,
                                    sigma_a, sigma_b,sigma_eps, sigma_eta, level)



