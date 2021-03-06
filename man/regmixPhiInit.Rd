% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regmix_funcs.R
\name{regmixPhiInit}
\alias{regmixPhiInit}
\title{regmixPhiInit}
\usage{
regmixPhiInit(y, x, z = NULL, parlist, h, tau, ninits = 1)
}
\arguments{
\item{y}{n by 1 vector of data for y}

\item{x}{n by q matrix of data for x}

\item{z}{n by p matrix of regressor associated with gam}

\item{parlist}{The parameter estimates as a list containing alpha, mu, sigma, and gam
in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
sigma = (sigma_1, ..., sigma_m), gam = (gam_1, ..., gam_m))}

\item{h}{h used as index for pivoting}

\item{tau}{Tau used to split the h-th component}

\item{ninits}{number of initial values to be generated}
}
\value{
A list with the following items:
\item{alpha}{m+1 by ninits matrix for alpha}
\item{mubeta}{q+1 by m+1 by ninits array for mu and beta}
\item{sigma}{m+1 by ninits matrix for sigma}
\item{gam}{m+1 by ninits matrix for gam}
}
\description{
Generates lists of parameters for initial candidates used by
the modified EM test for mixture of normals.
}

