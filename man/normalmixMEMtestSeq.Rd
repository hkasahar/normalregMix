\name{normalmixMEMtestSeq}
\alias{normalmixMEMtestSeq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Title %%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
normalmixMEMtestSeq(y, z = NULL, maxm = 3, ninits = 25, maxit = 2000,
                    nbtsp = 199, parallel = TRUE, cl = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{z}{
%%     ~~Describe \code{z} here~~
}
  \item{maxm}{
%%     ~~Describe \code{maxm} here~~
}
  \item{ninits}{
%%     ~~Describe \code{ninits} here~~
}
  \item{maxit}{
%%     ~~Describe \code{maxit} here~~
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
function (y, z = NULL, maxm = 3, ninits = 25, maxit = 2000, nbtsp = 199, 
    parallel = TRUE, cl = NULL) 
{
    y <- as.vector(y)
    n <- length(y)
    p <- 0
    if (!is.null(z)) {
        z <- as.matrix(z)
        p <- ncol(z)
        gamma <- matrix(0, nrow = p, ncol = maxm)
    }
    out <- vector("list", length = maxm)
    aic <- bic <- double(maxm)
    pvals <- emstat <- matrix(0, nrow = maxm, ncol = 3)
    loglik <- penloglik <- double(maxm)
    alpha <- mu <- sigma <- matrix(0, nrow = maxm, ncol = maxm)
    for (m in 1:maxm) {
        par0 <- normalmixPMLE(y = y, m = m, z = z, vcov.method = "none", 
            ninits = ninits, maxit = maxit)
        loglik[m] <- loglik0 <- par0$loglik
        penloglik[m] <- penloglik0 <- par0$penloglik
        aic[m] <- par0$aic
        bic[m] <- par0$bic
        parlist <- par0$parlist
        alpha0 <- parlist$alpha
        mu0 <- parlist$mu
        sigma0 <- parlist$sigma
        gamma0 <- parlist$gamma
        alpha[, m] <- c(alpha0, double(maxm - m))
        mu[, m] <- c(mu0, double(maxm - m))
        sigma[, m] <- c(sigma0, double(maxm - m))
        cat(sprintf("\%d-component model estimate:\n", m))
        tab = as.matrix(rbind(alpha0, mu0, sigma0))
        rownames(tab) <- c("alpha", "mu", "sigma")
        colnames(tab) <- c(paste("comp", ".", 1:m, sep = ""))
        print(tab, digits = 4)
        if (!is.null(z)) {
            gamma[, m] <- gamma0
            cat("gamma =", gamma0, "\n")
        }
        cat(sprintf("\nAIC, BIC, and loglike of 1 to \%.i", m), 
            "component models \n")
        cat(c("AIC    =", sprintf(" \%.2f", aic[1:m])), "\n")
        cat(c("BIC    =", sprintf(" \%.2f", bic[1:m])), "\n")
        cat(c("loglik =", sprintf("\%.2f", loglik[1:m])), "\n\n")
        cat(sprintf("Testing the null hypothesis of \%d components\n", 
            m))
        an <- anFormula(par = parlist, m = m, n = n)
        par1 <- normalmixMaxPhi(y = y, par0 = parlist, z = z, 
            an = an, ninits = ninits, maxit = maxit)
        emstat.m <- 2 * (par1$loglik - loglik0)
        cat(c("modified EM-test statitic ", sprintf("\%.3f", emstat.m)), 
            "\n")
        if (m <= 3) {
            em.out <- normalmixCrit(y = y, parlist = parlist, 
                z = z, values = emstat.m)
            cat(c("asymptotic p-values       ", sprintf("\%.3f", 
                em.out$pvals)), "\n \n")
        }
        else {
            em.out <- normalmixCritBoot(y = y, parlist = parlist, 
                z = z, values = emstat.m, ninits = ninits, nbtsp = nbtsp, 
                parallel = parallel, cl = cl)
            cat(c("bootstrap p-values        ", sprintf("\%.3f", 
                em.out$pvals)), "\n \n")
        }
        pvals[m, ] <- em.out$pvals
        emstat[m, ] <- emstat.m
    }
    a = list(alpha = alpha, mu = mu, sigma = sigma, gamma = gamma, 
        emstat = emstat, pvals = pvals, aic = aic, bic = bic, 
        loglik = loglik, penloglik = penloglik)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
