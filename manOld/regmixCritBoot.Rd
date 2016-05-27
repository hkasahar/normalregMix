\name{regmixCritBoot}
\alias{regmixCritBoot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
title %%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
regmixCritBoot(y, x, parlist, z = NULL, values = NULL, ninits = 100,
                nbtsp = 199, parallel = TRUE, cl = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{parlist}{
%%     ~~Describe \code{parlist} here~~
}
  \item{z}{
%%     ~~Describe \code{z} here~~
}
  \item{values}{
%%     ~~Describe \code{values} here~~
}
  \item{ninits}{
%%     ~~Describe \code{ninits} here~~
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
function (y, x, parlist, z = NULL, values = NULL, ninits = 100, 
    nbtsp = 199, parallel = TRUE, cl = NULL) 
{
    set.seed(123456)
    y <- as.vector(y)
    n <- length(y)
    x <- as.matrix(x)
    alpha <- parlist$alpha
    mubeta <- parlist$mubeta
    sigma <- parlist$sigma
    gamma <- parlist$gamma
    m <- length(alpha)
    pvals <- NULL
    ybset <- replicate(nbtsp, rnormregmix(n = n, x = x, alpha = alpha, 
        mubeta = mubeta, sigma = sigma))
    if (!is.null(z)) {
        zgamma <- as.matrix(z) \%*\% gamma
        ybset <- ybset + as.vector(zgamma)
    }
    if (parallel) {
        if (is.null(cl)) {
            ncpus <- parallel::detectCores()
            cl <- parallel::makePSOCKcluster(rep("localhost", 
                ncpus))
            parallel::clusterSetRNGStream(cl, 123456)
            out <- parallel::parLapplyLB(cl, 1:nbtsp, function(j) regmixMEMtest(y = ybset[, 
                j], x = x, m = m, z = z, ninits = ninits))
            parallel::stopCluster(cl)
        }
        else {
            parallel::clusterSetRNGStream(cl, 123456)
            out <- parallel::parLapplyLB(cl, 1:nbtsp, function(j) regmixMEMtest(y = ybset[, 
                j], x = x, m = m, z = z, ninits = ninits))
        }
    }
    else {
        out <- apply(ybset, 2, regmixMEMtest, x = x, m = m, z = z, 
            ninits = ninits)
    }
    emstat.b <- sapply(out, "[[", "emstat")
    emstat.b <- t(apply(emstat.b, 1, sort))
    q <- ceiling(nbtsp * c(0.9, 0.95, 0.99))
    crit <- emstat.b[, q]
    if (!is.null(values)) {
        pvals <- rowMeans(emstat.b > values)
    }
    return(list(crit = crit, pvals = pvals))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
