\name{normalregMix-package}
\alias{normalregMix-package}
\alias{normalregMix}
\docType{package}
\title{
modified EM test for the number of components in normal mixtures
and normal regression mixtures
}
\description{
More about what it does (maybe more than one line)
~~ A concise (1-5 lines) description of the package ~~
}
\details{
\tabular{ll}{
Package: \tab normalregMix\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-04-07\cr
License: \tab What license is it under?\cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Hiroyuki Kasahara and Katsumi Shimotsu
Maintainer: Katsumi Shimotsu <shimotsu@e.u-tokyo.ac.jp>
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
alpha <- c(0.5, 0.5)
mu    <- c(-1, 1)
sigma <- c(1, 1)
y <- rnormregmix(100, , alpha, mu, sigma)
normalmixPMLE(y, m = 2)
normalmixMEMtest(y, m = 2)
}
