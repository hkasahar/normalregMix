\name{regmixPMLE}
\alias{regmixPMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Penalized Maximum Likelihood Estimation of a Normal Mixture Regression Model
}
\description{
Estimate parameters of a normal mixture regression model by the method of penalized maximum likelhood.
}
\usage{
regmixPMLE(y, x, m = 2, z = NULL, vcov.method = c("Hessian", "OPG", "none"),
            ninits = 100, epsilon = 1e-08, maxit = 2000, epsilon.short = 0.01,
            maxit.short = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Numeric vector of dependent variable.
}
  \item{x}{
An nxq matrix of explanatory variables. Should not include an intercept.
}
  \item{m}{
The number of components in the mixture.}
  \item{z}{
An nxp matrix of covariates whose coefficients are common to all the components.
}
  \item{vcov.method}{
The method to compute the covariance matrix of the PMLE
}
  \item{ninits}{
The number of randomly drawn initial values.
}
  \item{epsilon}{
The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}. 
}
  \item{maxit}{
The maximum number of iterations.
}
  \item{epsilon.short}{
The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
}
  \item{maxit.short}{
The maximum number of iterations in short EM.
}
}
\details{
The likelihood function of normal mixtures diverges to infinity when the variance of a component is allowed to go to zero. The penalized MLE overcomes this problem by the use of a penalty term. See Chen et al. (2008) for details.
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
function (y, x, m = 2, z = NULL, vcov.method = c("Hessian", "OPG", 
    "none"), ninits = 100, epsilon = 1e-08, maxit = 2000, epsilon.short = 0.01, 
    maxit.short = 500) 
{
    y <- as.vector(y)
    x <- as.matrix(x)
    n <- length(y)
    if (nrow(x) != n) {
        stop("y and x must have the same number of rows.")
    }
    x1 <- cbind(1, x)
    k <- ncol(x1)
    p <- 0
    gamma <- NULL
    ninits.short <- ninits * 4 * (k + p) * m
    vcov.method <- match.arg(vcov.method)
    if (!is.null(z)) {
        z <- as.matrix(z)
        p <- ncol(z)
        if (nrow(z) != n) {
            stop("y and z must have the same number of rows.")
        }
    }
    npar <- m - 1 + (k + 1) * m + p
    xz <- cbind(x, z)
    ls.out <- lsfit(xz, y)
    sd0 <- sqrt(mean(ls.out$residuals^2))
    if (m == 1) {
        mubeta <- as.matrix(unname(ls.out$coeff[1:k]))
        if (!is.null(z)) {
            gamma <- unname(ls.out$coeff[(k + 1):(k + p)])
        }
        res <- ls.out$residuals
        sigma <- sqrt(mean(res * res))
        loglik <- -(n/2) * (1 + log(2 * pi) + 2 * log(sigma))
        aic <- -2 * loglik + 2 * npar
        bic <- -2 * loglik + log(n) * npar
        penloglik <- loglik
        parlist <- list(alpha = 1, mubeta = mubeta, sigma = sigma, 
            gamma = gamma)
        coefficients <- c(alpha = 1, mubeta = mubeta, sigma = sigma, 
            gamma = gamma)
        postprobs <- rep(1, n)
    }
    else {
        tmp <- regmixPMLEinit(y = y, x = x, z = z, ninits = ninits.short, 
            m = m)
        alphaset.s <- tmp$alpha
        mubetaset.s <- matrix(tmp$mubeta, nrow = m * k)
        sigmaset.s <- tmp$sigma
        gammaset.s <- tmp$gamma
        h <- 0
        sigma0 <- rep(sd0, m)
        mu0 <- double(m)
        an <- 1/n
        jpvt <- integer(max(k, p))
        if (is.null(z)) {
            setting <- c(n, m, k, ninits.short, maxit.short, 
                jpvt)
            out.short <- .C("regmixpmle", as.integer(setting), 
                as.double(y), as.double(x1), alphaset = as.double(alphaset.s), 
                mubetaset = as.double(mubetaset.s), sigmaset = as.double(sigmaset.s), 
                as.double(sigma0), as.double(mu0), as.double(an), 
                as.integer(h), lub = double(2 * m), double(3 * 
                  m), post = double(n * m), loglikset = double(ninits.short), 
                penloglikset = double(ninits.short), notcg = integer(ninits.short), 
                as.double(epsilon.short), double(n * (k + 1)), 
                package = "normalregMix")
        }
        else {
            setting.z <- c(n, m, k, p, ninits.short, maxit.short, 
                jpvt)
            out.short <- .C("regmixpmle_z", as.integer(setting.z), 
                as.double(y), as.double(x1), as.double(z), alphaset = as.double(alphaset.s), 
                mubetaset = as.double(mubetaset.s), sigmaset = as.double(sigmaset.s), 
                gammaset = as.double(gammaset.s), as.double(sigma0), 
                as.double(mu0), as.double(an), as.integer(h), 
                lub = double(2 * m), double(3 * m), post = double(n * 
                  m), loglikset = double(ninits.short), penloglikset = double(ninits.short), 
                notcg = integer(ninits.short), as.double(epsilon.short), 
                double(n * (k + 3 + p)), package = "normalregMix")
        }
        penloglik.short <- out.short$penloglikset
        oo <- order(penloglik.short, decreasing = TRUE)
        oo.inits <- oo[1:ninits]
        alphaset <- alphaset.s[, oo.inits]
        mubetaset <- mubetaset.s[, oo.inits]
        sigmaset <- sigmaset.s[, oo.inits]
        if (is.null(z)) {
            setting <- c(n, m, k, ninits, maxit, jpvt)
            out <- .C("regmixpmle", as.integer(setting), as.double(y), 
                as.double(x1), alphaset = as.double(alphaset), 
                mubetaset = as.double(mubetaset), sigmaset = as.double(sigmaset), 
                as.double(sigma0), as.double(mu0), as.double(an), 
                as.integer(h), lub = double(2 * m), double(3 * 
                  m), post = double(n * m), loglikset = double(ninits), 
                penloglikset = double(ninits), notcg = integer(ninits), 
                as.double(epsilon), double(n * (k + 1)), package = "normalregMix")
        }
        else {
            setting.z <- c(n, m, k, p, ninits, maxit, jpvt)
            gammaset <- gammaset.s[, oo.inits]
            out <- .C("regmixpmle_z", as.integer(setting.z), 
                as.double(y), as.double(x1), as.double(z), alphaset = as.double(alphaset), 
                mubetaset = as.double(mubetaset), sigmaset = as.double(sigmaset), 
                gammaset = as.double(gammaset), as.double(sigma0), 
                as.double(mu0), as.double(an), as.integer(h), 
                lub = double(2 * m), double(3 * m), post = double(n * 
                  m), loglikset = double(ninits), penloglikset = double(ninits), 
                notcg = integer(ninits), as.double(epsilon), 
                double(n * (k + 3 + p)), package = "normalregMix")
        }
        if (mean(out$notcg) >= 0.9) {
            warning(sprintf("The EM algorithm failed to converge in \%d\%\% of the initial values.\n Try increasing maxit.", 
                100 * mean(out$notcg)))
        }
        index <- which.max(out$penloglikset)
        alpha <- out$alphaset[(m * (index - 1) + 1):(m * index)]
        mubeta <- matrix(out$mubetaset[(k * m * (index - 1) + 
            1):(k * m * index)], nrow = k)
        sigma <- out$sigmaset[(m * (index - 1) + 1):(m * index)]
        gamma <- out$gammaset[(p * (index - 1) + 1):(p * index)]
        penloglik <- out$penloglikset[index]
        loglik <- out$loglikset[index]
        if (is.null(z)) {
            setting <- c(n, m, k, 1, 1, jpvt)
            out.p <- .C("regmixpmle", as.integer(setting), as.double(y), 
                as.double(x1), alphaset = as.double(alpha), mubetaset = as.double(mubeta), 
                sigmaset = as.double(sigma), as.double(sigma0), 
                as.double(mu0), as.double(an), as.integer(h), 
                lub = double(2 * m), double(3 * m), post = double(n * 
                  m), loglikset = double(1), penloglikset = double(1), 
                notcg = integer(1), as.double(epsilon), double(n * 
                  (k + 1)), package = "normalregMix")
        }
        else {
            setting.z <- c(n, m, k, p, 1, 1, jpvt)
            out.p <- .C("regmixpmle_z", as.integer(setting.z), 
                as.double(y), as.double(x1), as.double(z), alphaset = as.double(alpha), 
                mubetaset = as.double(mubeta), sigmaset = as.double(sigma), 
                gammaset = as.double(gamma), as.double(sigma0), 
                as.double(mu0), as.double(an), as.integer(h), 
                lub = double(2 * m), double(3 * m), post = double(n * 
                  m), loglikset = double(1), penloglikset = double(1), 
                notcg = integer(1), as.double(epsilon), double(n * 
                  (k + 3 + p)), package = "normalregMix")
        }
        aic <- -2 * loglik + 2 * npar
        bic <- -2 * loglik + log(n) * npar
        mu.order <- order(mubeta[1, ])
        alpha <- alpha[mu.order]
        mubeta <- mubeta[, mu.order]
        sigma <- sigma[mu.order]
        mu <- mubeta[1, ]
        beta <- mubeta[-1, ]
        postprobs <- matrix(out.p$post, nrow = n)
        postprobs <- postprobs[, mu.order]
        colnames(postprobs) <- c(paste("comp", ".", 1:m, sep = ""))
        mubeta.name <- matrix(0, nrow = k, ncol = m)
        mubeta.name[1, ] <- paste("mu", 1:m, sep = "")
        for (i in 1:(k - 1)) {
            for (j in 1:m) {
                mubeta.name[i + 1, j] <- paste("beta", j, i, 
                  sep = "")
            }
        }
        parlist <- list(alpha = alpha, mubeta = mubeta, sigma = sigma, 
            gamma = gamma)
        coefficients <- unlist(parlist)
        names(coefficients)[(m + 1):((k + 1) * m)] <- c(mubeta.name)
    }
    if (vcov.method == "none") {
        vcov <- NULL
    }
    else {
        vcov <- regmixVcov(y = y, x = x, coefficients = coefficients, 
            z = z, vcov.method = vcov.method)
    }
    a <- list(coefficients = coefficients, parlist = parlist, 
        vcov = vcov, loglik = loglik, penloglik = penloglik, 
        aic = aic, bic = bic, postprobs = postprobs, call = match.call(), 
        m = m)
    class(a) <- "normalregMix"
    a
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
