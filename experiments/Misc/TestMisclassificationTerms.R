# Returns a misclassification rate given two components i, j
# where omega_JI denotes the probability of classifying as jth component
# where the true model is ith component 
GetOmegaJI <- function(phi_i, phi_j) {
  alpha_i <- phi_i[[1]]
  alpha_j <- phi_j[[1]]
  mu_i <- phi_i[[2]]
  mu_j <- phi_j[[2]]
  sigma_i <- phi_i[[3]]
  sigma_j <- phi_j[[3]]
  
  a <- (1/sigma_j^2 - 1/sigma_i^2)
  b <- mu_i / sigma_i^2 - mu_j / sigma_j^2
  c <- mu_j^2 / sigma_j^2 - mu_i^2 / sigma_i^2

  # WARNING: What is normp? according to one MATLAB package description,
  # it is the cdf of a standard normal. If that's the case, this should work.  
  if (sigma_i == sigma_j)
    if (mu_i > mu_j)
      omega_ji = pnorm((2 * log(alpha_j/alpha_i) - c)/(2*b), 
                       mean = mu_i, sd = sigma_i)
    else
      omega_ji = 1 - pnorm((2 * log(alpha_j/alpha_i) - c)/(2*b), 
                           mean = mu_i, sd = sigma_i)
  else {
    d <- 2 * log(alpha_j * sigma_i / (alpha_i * sigma_j)) - c + (b^2 / a)
    da <- max(d/a, 0)
    if (sigma_i > sigma_j)
      omega_ji = pnorm(sqrt(da)-b/a, mean = mu_i, sd = sigma_i) -
                 pnorm(-sqrt(da)-b/a, mean = mu_i, sd = sigma_i)
    else 
      omega_ji = 1 + 
                  pnorm(-sqrt(da)-b/a, mean = mu_i, sd = sigma_i) -
                  pnorm(sqrt(da)-b/a, mean = mu_i, sd = sigma_i)
  }
  return (omega_ji)
}
GetSimulatedOmegaJI <- function(phi_i, phi_j) {
  alpha_i <- phi_i[[1]]
  alpha_j <- phi_j[[1]]
  mu_i <- phi_i[[2]]
  mu_j <- phi_j[[2]]
  sigma_i <- phi_i[[3]]
  sigma_j <- phi_j[[3]]
  sample <- rnorm(30000, mean = mu_i, sd = sigma_i)
  mean(alpha_i * dnorm(sample, mean = mu_i, sd = sigma_i) < 
         alpha_j * dnorm(sample, mean = mu_j, sd = sigma_j))
}
# Case 1 (sigma_i = sigma_j, mu_i > mu_j)
phi_i <- list(alpha = 0.6, mu = 3, sigma = 2)
phi_j <- list(alpha = 0.4, mu = 1, sigma = 2)
c(GetSimulatedOmegaJI(phi_i, phi_j), GetOmegaJI(phi_i, phi_j))

# Case 2 (sigma_i = sigma_j, mu_i < mu_j)
phi_i <- list(alpha = 0.6, mu = 1, sigma = 2)
phi_j <- list(alpha = 0.4, mu = 3, sigma = 2)
c(GetSimulatedOmegaJI(phi_i, phi_j), GetOmegaJI(phi_i, phi_j))

# Case 3 (sigma_i > sigma_j)
phi_i <- list(alpha = 0.6, mu = 3, sigma = 4)
phi_j <- list(alpha = 0.4, mu = 1, sigma = 2)
c(GetSimulatedOmegaJI(phi_i, phi_j), GetOmegaJI(phi_i, phi_j))
phi_i <- list(alpha = 0.6, mu = 1, sigma = 4)
phi_j <- list(alpha = 0.4, mu = 3, sigma = 2)
c(GetSimulatedOmegaJI(phi_i, phi_j), GetOmegaJI(phi_i, phi_j))

# Case 4 (sigma_i < sigma_j)
phi_i <- list(alpha = 0.6, mu = 3, sigma = 1)
phi_j <- list(alpha = 0.4, mu = 1, sigma = 5)
c(GetSimulatedOmegaJI(phi_i, phi_j), GetOmegaJI(phi_i, phi_j))
phi_i <- list(alpha = 0.6, mu = 1, sigma = 1)
phi_j <- list(alpha = 0.4, mu = 3, sigma = 5)
c(GetSimulatedOmegaJI(phi_i, phi_j), GetOmegaJI(phi_i, phi_j))
