\name{regmixCrit}
\alias{regmixCrit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
title %%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
regmixCrit(y, x, parlist, z = NULL, values = NULL, parallel = TRUE,
            cl = NULL, nrep = 1000, ninits.crit = 25)
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
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
}
  \item{cl}{
%%     ~~Describe \code{cl} here~~
}
  \item{nrep}{
%%     ~~Describe \code{nrep} here~~
}
  \item{ninits.crit}{
%%     ~~Describe \code{ninits.crit} here~~
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
function (y, x, parlist, z = NULL, values = NULL, parallel = TRUE, 
    cl = NULL, nrep = 1000, ninits.crit = 25) 
{
    y <- as.vector(y)
    n <- length(y)
    p <- 0
    alpha <- parlist$alpha
    mubeta <- parlist$mubeta
    sigma <- parlist$sigma
    gamma <- parlist$gamma
    m <- length(alpha)
    if (!is.null(z)) {
        z <- as.matrix(z)
        p <- ncol(z)
        y <- y - z \%*\% gamma
    }
    pvals <- NULL
    x <- as.matrix(x)
    x1 <- cbind(1, x)
    q <- ncol(x)
    W <- t(t(y - x1 \%*\% mubeta)/sigma)
    f <- t(t(exp(-W^2/2)/sqrt(2 * pi))/sigma)
    f0 <- colSums(t(f) * alpha)
    H <- hermite(W, sigma)
    w_m <- t(t(f) * alpha)/f0
    if (m == 1) {
        S_alpha <- NULL
    }
    else {
        S_alpha <- (f[, 1:(m - 1)] - f[, m])/f0
    }
    S_mu <- w_m * H[, , 1]
    S_beta <- matrix(0, nrow = n, ncol = q * m)
    for (j in 1:m) {
        S_beta[, (1 + (j - 1) * q):(j * q)] <- x * S_mu[, j]
    }
    S_sigma <- w_m * H[, , 2]
    S_eta <- cbind(S_alpha, S_mu, S_beta, S_sigma)
    if (!is.null(z)) {
        S_gamma <- z * rowSums(S_mu)
        S_eta <- cbind(S_eta, S_gamma)
    }
    n_lam <- q * (q + 1)/2 + 2 * q + 2
    S_lam <- matrix(0, nrow = n, ncol = n_lam * m)
    xx <- matrix(0, nrow = n, ncol = q * (q + 1)/2)
    xx[, 1:q] <- x * x
    if (q > 1) {
        t <- q + 1
        for (j in 1:(q - 1)) {
            for (i in (j + 1):q) {
                xx[, t] <- 2 * x[, j] * x[, i]
                t <- t + 1
            }
        }
    }
    for (j in 1:m) {
        w_2 <- S_sigma[, j]
        w_3 <- w_m[, j] * H[, j, 3]
        w_4 <- w_m[, j] * H[, j, 4]
        S_lam_1 <- cbind(w_3, x * w_2)
        S_lam_2 <- cbind(w_4, x * w_3, xx * w_2)
        S_lam[, ((j - 1) * n_lam + 1):(j * n_lam)] <- cbind(S_lam_1, 
            S_lam_2)
    }
    I_eta <- t(S_eta) \%*\% S_eta/n
    I_lam <- t(S_lam) \%*\% S_lam/n
    I_el <- t(S_eta) \%*\% S_lam/n
    if (qr(I_eta)$rank == nrow(I_eta)) {
        I_eta_lam <- I_lam - t(I_el) \%*\% solve(I_eta, I_el)
    }
    else {
        stop("The critical value cannot be computed due to singularity of some matrices.\n     Please try a bootstrap version, regmixCritBoot and regmixMEMtestBoot.")
    }
    set.seed(123456)
    e <- eigen(I_eta_lam, symmetric = TRUE)
    u <- t(e$vec \%*\% (t(e$vec) * sqrt(e$val)) \%*\% matrix(rnorm(nrep * 
        n_lam * m), nrow = n_lam * m))
    q_1 <- 1 + q
    q_2 <- 1 + q + q * (q + 1)/2
    LR <- matrix(0, nrow = nrep, ncol = m)
    if ((parallel) && (is.null(cl))) {
        ncpus <- parallel::detectCores()
        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
    }
    for (j in 1:m) {
        I_j <- I_eta_lam[((j - 1) * n_lam + 1):(j * n_lam), ((j - 
            1) * n_lam + 1):(j * n_lam)]
        if (qr(I_j)$rank == nrow(I_j)) {
            Z_j <- u[, ((j - 1) * n_lam + 1):(j * n_lam)] \%*\% 
                solve(I_j)
        }
        else {
            stop("The critical value cannot be computed due to singularity of some matrices.\n     Please try a bootstrap version, regmixCritBoot and regmixMEMtestBoot.")
        }
        V <- solve(I_j)
        n_v <- ncol(V)
        V_11 <- V[1:q_1, 1:q_1]
        V_12 <- V[1:q_1, (q_1 + 1):n_v]
        V_21 <- t(V_12)
        V_22 <- V[(q_1 + 1):n_v, (q_1 + 1):n_v]
        Z_1 <- Z_j[, 1:q_1]
        Z_2 <- Z_j[, (q_1 + 1):ncol(Z_j)]
        V_1_2 <- V_11 - V_12 \%*\% solve(V_22, V_21)
        Z_1_2 <- Z_1 - Z_2 \%*\% solve(V_22, V_21)
        inv_V_22 <- solve(V_22)
        Z_22 <- t(inv_V_22[1, ] \%*\% t(Z_2))
        LR_1 <- rowSums((Z_1_2 \%*\% solve(V_1_2)) * Z_1_2) + (1/inv_V_22[1, 
            1]) * (Z_22^2) * (Z_22 < 0)
        if (parallel) {
            parallel::clusterSetRNGStream(cl, 123456)
            LR_2 <- parallel::parRapply(cl, Z_j, LR_2.comp, I_j, 
                q, ninits.crit)
        }
        else {
            LR_2 <- apply(Z_j, 1, LR_2.comp, I_j, q, ninits.crit)
        }
        LR[, j] <- apply(cbind(LR_1, LR_2), 1, max)
    }
    if (parallel) {
        parallel::stopCluster(cl)
    }
    max_EM <- apply(LR, 1, max)
    max_EM_sort <- sort(max_EM)
    qc <- floor(nrep * c(0.9, 0.95, 0.99))
    crit <- max_EM_sort[qc]
    if (!is.null(values)) {
        k <- length(values)
        pvals <- rowMeans(t(matrix(rep.int(max_EM_sort, k), ncol = k)) >= 
            values)
    }
    return(list(crit = crit, pvals = pvals))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
