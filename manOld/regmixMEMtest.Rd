\name{regmixMEMtest}
\alias{regmixMEMtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
title%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
regmixMEMtest(y, x, m = 2, z = NULL, an = NULL, ninits = 100,
              crit.method = c("none", "asy", "boot"), parallel = TRUE,
              cl = NULL, nbtsp = 199)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{z}{
%%     ~~Describe \code{z} here~~
}
  \item{an}{
%%     ~~Describe \code{an} here~~
}
  \item{ninits}{
%%     ~~Describe \code{ninits} here~~
}
  \item{crit.method}{
%%     ~~Describe \code{crit.method} here~~
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
}
  \item{cl}{
%%     ~~Describe \code{cl} here~~
}
  \item{nbtsp}{
%%     ~~Describe \code{nbtsp} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, x, m = 2, z = NULL, an = NULL, ninits = 100, crit.method = c("none", 
    "asy", "boot"), parallel = TRUE, cl = NULL, nbtsp = 199) 
{
    y <- as.vector(y)
    x <- as.matrix(x)
    q <- ncol(x)
    crit.method <- match.arg(crit.method)
    if (is.null(an)) {
        an <- switch(as.character(q), `1` = 2.2, `2` = 3.6, `3` = 4.6, 
            `4` = 8.2, 2.2)
    }
    if (!is.null(z)) {
        z <- as.matrix(z)
    }
    par0 <- regmixPMLE(y = y, x = x, m = m, z = z, vcov.method = "none", 
        ninits = ninits)
    loglik0 <- par0$loglik
    par1 <- regmixMaxPhi(y = y, x = x, par0 = par0$parlist, z = z, 
        an = an, ninits = ninits)
    emstat <- 2 * (par1$loglik - loglik0)
    if (crit.method == "asy") {
        result <- regmixCrit(y = y, x = x, parlist = par0$parlist, 
            z = z, values = emstat, parallel = parallel, cl = cl, 
            nrep = 1000, ninits.crit = 25)
    }
    else if (crit.method == "boot") {
        result <- regmixCritBoot(y = y, x = x, parlist = par0$parlist, 
            z = z, values = emstat, ninits = ninits, nbtsp = nbtsp, 
            parallel = parallel, cl = cl)
    }
    else {
        result <- list()
        result$crit <- result$pvals <- rep(NA, 3)
    }
    a <- list(emstat = emstat, pvals = result$pvals, crit = result$crit, 
        parlist = par0$parlist, call = match.call(), m = m, crit.method = crit.method, 
        nbtsp = nbtsp)
    class(a) <- "normalregMix"
    a
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
