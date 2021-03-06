% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalmix_funcs.R
\name{normalmixPMLEinit}
\alias{normalmixPMLEinit}
\title{normalmixPMLEinit}
\usage{
normalmixPMLEinit(y, z = NULL, ninits = 1, m = 2)
}
\arguments{
\item{y}{n by 1 vector of data}

\item{z}{n by p matrix of regressor associated with gamma}

\item{ninits}{number of initial values to be generated}

\item{m}{The number of components in the mixture}
}
\value{
A list with the following items:
\item{alpha}{m by ninits matrix for alpha}
\item{mu}{m by ninits matrix for mu}
\item{sigma}{m by ninits matrix for sigma}
\item{gam}{m by ninits matrix for gam}
}
\description{
Generate initial values used by the PMLE of mixture of normals
}

