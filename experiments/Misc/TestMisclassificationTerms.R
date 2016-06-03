GetSimulatedOmegaJI <- function(phi_i, phi_j) {
  alpha_i <- phi_i[1]
  alpha_j <- phi_j[1]
  mu_i <- phi_i[2]
  mu_j <- phi_j[2]
  sigma_i <- phi_i[3]
  sigma_j <- phi_j[3]
  sample <- rnorm(30000, mean = mu_i, sd = sigma_i)
  mean(alpha_i * dnorm(sample, mean = mu_i, sd = sigma_i) < 
         alpha_j * dnorm(sample, mean = mu_j, sd = sigma_j))
}
# Case 1 (sigma_i = sigma_j, mu_i > mu_j)
phi_i <- c(alpha = 0.6, mu = 3, sigma = 2)
phi_j <- c(alpha = 0.4, mu = 1, sigma = 2)
c(simulated = GetSimulatedOmegaJI(phi_i, phi_j), actual = omega.ji(phi_i, phi_j))

# Case 2 (sigma_i = sigma_j, mu_i < mu_j)
phi_i <- c(alpha = 0.6, mu = 1, sigma = 2)
phi_j <- c(alpha = 0.4, mu = 3, sigma = 2)
c(simulated = GetSimulatedOmegaJI(phi_i, phi_j), actual = omega.ji(phi_i, phi_j))

# Case 3 (sigma_i > sigma_j)
phi_i <- c(alpha = 0.6, mu = 3, sigma = 4)
phi_j <- c(alpha = 0.4, mu = 1, sigma = 2)
c(simulated = GetSimulatedOmegaJI(phi_i, phi_j), actual = omega.ji(phi_i, phi_j))
phi_i <- c(alpha = 0.6, mu = 1, sigma = 4)
phi_j <- c(alpha = 0.4, mu = 3, sigma = 2)
c(simulated = GetSimulatedOmegaJI(phi_i, phi_j), actual = omega.ji(phi_i, phi_j))

# Case 4 (sigma_i < sigma_j)
phi_i <- c(alpha = 0.6, mu = 3, sigma = 1)
phi_j <- c(alpha = 0.4, mu = 1, sigma = 5)
c(simulated = GetSimulatedOmegaJI(phi_i, phi_j), actual = omega.ji(phi_i, phi_j))
phi_i <- c(alpha = 0.6, mu = 1, sigma = 1)
phi_j <- c(alpha = 0.4, mu = 3, sigma = 5)
c(simulated = GetSimulatedOmegaJI(phi_i, phi_j), actual = omega.ji(phi_i, phi_j))
