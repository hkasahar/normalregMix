\name{normalmixCrit}
\alias{normalmixCrit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes critical values and p-values of the modified EM test of normal mixtures%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
normalmixCrit(y, parlist, z = NULL, values = NULL, nrep = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
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
  \item{nrep}{
%%     ~~Describe \code{nrep} here~~
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
function (y, parlist, z = NULL, values = NULL, nrep = 10000) 
{
    y <- as.vector(y)
    n <- length(y)
    p <- 0
    alpha <- parlist$alpha
    mu <- parlist$mu
    sigma <- parlist$sigma
    gamma <- parlist$gamma
    m <- length(alpha)
    if (!is.null(z)) {
        z <- as.matrix(z)
        p <- ncol(z)
        y <- y - z \%*\% gamma
    }
    pvals <- NULL
    if (m == 1) {
        crit <- qchisq(c(0.1, 0.05, 0.01), 2, lower.tail = F)
        if (!is.null(values)) {
            k <- length(values)
            pvals <- pchisq(values, 2, lower.tail = F)
        }
    }
    else {
        set.seed(123456)
        Z0 <- t((t(matrix(rep.int(y, m), ncol = m)) - mu)/sigma)
        f <- t(t(exp(-Z0^2/2)/sqrt(2 * pi))/sigma)
        f0 <- colSums(t(f) * alpha)
        S_alpha <- (f[, 1:(m - 1)] - f[, m])/f0
        H <- hermite(Z0, sigma)
        S_mu <- t(t(H[, , 1] * f) * alpha)/f0
        S_sigma <- t(t(H[, , 2] * f) * alpha)/f0
        S_gamma <- rowSums(S_mu) * z
        S_lambda21 <- t(t(H[, , 3] * f) * alpha)/f0
        S_lambda22 <- t(t(H[, , 4] * f) * alpha)/f0
        S_lambda <- matrix(rbind(S_lambda21, S_lambda22), nrow = n, 
            ncol = 2 * m)
        S_eta <- cbind(S_alpha, S_gamma, S_mu, S_sigma)
        I11 <- t(S_eta) \%*\% S_eta/n
        I21 <- t(S_lambda) \%*\% S_eta/n
        I22 <- t(S_lambda) \%*\% S_lambda/n
        Iall <- rbind(cbind(I11, t(I21)), cbind(I21, I22))
        if (rcond(Iall) < .Machine$double.eps) {
            eig <- eigen(Iall, symmetric = TRUE)
            tol2 <- (1e-14) * eig$values[1]
            vals <- eig$values
            vals[vals < tol2] <- tol2
            Iall.mod <- eig$vectors \%*\% (vals * t(eig$vectors))
            I11 <- Iall.mod[1:(p + 3 * m - 1), 1:(p + 3 * m - 
                1)]
            I21 <- Iall.mod[(p + 3 * m):(p + 5 * m - 1), 1:(p + 
                3 * m - 1)]
            I22 <- Iall.mod[(p + 3 * m):(p + 5 * m - 1), (p + 
                3 * m):(p + 5 * m - 1)]
        }
        I221 <- I22 - I21 \%*\% solve(I11, t(I21))
        e <- eigen(I221, symmetric = TRUE)
        u <- t(e$vec \%*\% (t(e$vec) * sqrt(e$val)) \%*\% matrix(rnorm(nrep * 
            2 * m), 2 * m, nrep))
        EM <- matrix(0, nrow = nrep, ncol = m)
        for (jj in 1:m) {
            uu <- u[, (2 * jj - 1):(2 * jj)]
            Ivinv <- solve(I221[(2 * jj - 1):(2 * jj), (2 * jj - 
                1):(2 * jj)])
            EM[, jj] <- rowSums((uu \%*\% Ivinv) * uu)
        }
        max_EM <- apply(EM, 1, max)
        max_EM_sort <- sort(max_EM)
        q <- ceiling(nrep * c(0.9, 0.95, 0.99))
        crit <- max_EM_sort[q]
        if (!is.null(values)) {
            k <- length(values)
            pvals <- rowMeans(t(matrix(rep.int(max_EM_sort, k), 
                ncol = k)) >= values)
        }
    }
    return(list(crit = crit, pvals = pvals))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
