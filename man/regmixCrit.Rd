% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regmix_crit.R
\name{regmixCrit}
\alias{regmixCrit}
\title{regmixCrit}
\usage{
regmixCrit(y, x, parlist, z = NULL, values = NULL, parallel = 1,
  cl = NULL, nrep = 1000, ninits.crit = 25)
}
\arguments{
\item{y}{n by 1 vector of data for y}

\item{x}{n by q matrix of data for x}

\item{parlist}{The parameter estimates as a list containing alpha, mu, sigma, and gamma
in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
sigma = (sigma_1, ..., sigma_m), gam = (gamma_1, ..., gamma_m))}

\item{z}{n by p matrix of regressor associated with gamma}

\item{values}{3 by 1 Vector of length 3 (k = 1, 2, 3) at which the p-values are computed}

\item{parallel}{Determines whether package \code{doParallel} is used for calculation}

\item{cl}{Cluster used for parallelization; if it is \code{NULL}, the system will automatically}

\item{nrep}{The number of replications used to compute p-values}

\item{ninits.crit}{The number of initial guesses to form critical values}
}
\value{
A list with the following items:
\item{crit}{3 by 3 matrix of (0.1, 0.05, 0.01 critical values), jth row corresponding to k=j}
\item{pvals}{A vector of p-values at k = 1, 2, 3}
}
\description{
Computes the critical values of the modified EM test.
}

