\name{normalmixMEMtest}
\alias{normalmixMEMtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
modified em test%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
normalmixMEMtest(y, m = 2, z = NULL, an = NULL, ninits = 25,
                  crit.method = c("asy", "boot", "none"),
                  nbtsp = 199, parallel = TRUE, cl = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
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
  \item{nbtsp}{
%%     ~~Describe \code{nbtsp} here~~
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
}
  \item{cl}{
%%     ~~Describe \code{cl} here~~
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
function (y, m = 2, z = NULL, an = NULL, ninits = 25, crit.method = c("asy", 
    "boot", "none"), nbtsp = 199, parallel = TRUE, cl = NULL) 
{
    y <- as.vector(y)
    n <- length(y)
    if (!is.null(z)) {
        z <- as.matrix(z)
    }
    crit.method <- match.arg(crit.method)
    par0 <- normalmixPMLE(y = y, m = m, z = z, vcov.method = "none", 
        ninits = ninits)
    loglik0 <- par0$loglik
    if (is.null(an)) {
        an <- anFormula(par = par0$parlist, m = m, n = n)
    }
    par1 <- normalmixMaxPhi(y = y, par0 = par0$parlist, z = z, 
        an = an, ninits = ninits)
    emstat <- 2 * (par1$loglik - loglik0)
    if (crit.method == "asy") {
        result <- normalmixCrit(y = y, parlist = par0$parlist, 
            z = z, values = emstat)
    }
    else if (crit.method == "boot") {
        result <- normalmixCritBoot(y = y, parlist = par0$parlist, 
            z = z, values = emstat, ninits = ninits, nbtsp = nbtsp, 
            parallel = parallel, cl = cl)
    }
    else {
        result <- list()
        result$crit <- result$pvals <- rep(NA, 3)
    }
    a <- list(emstat = emstat, pvals = result$pvals, crit = result$crit, 
        crit.method = crit.method, parlist = par0$parlist, call = match.call(), 
        m = m)
    class(a) <- "normalregMix"
    a
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
