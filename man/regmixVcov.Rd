\name{regmixVcov}
\alias{regmixVcov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
title %%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
regmixVcov(y, x, coefficients, z = NULL, vcov.method = c("Hessian", "OPG"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{coefficients}{
%%     ~~Describe \code{coefficients} here~~
}
  \item{z}{
%%     ~~Describe \code{z} here~~
}
  \item{vcov.method}{
%%     ~~Describe \code{vcov.method} here~~
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
function (y, x, coefficients, z = NULL, vcov.method = c("Hessian", 
    "OPG")) 
{
    y <- as.vector(y)
    n <- length(y)
    len <- length(coefficients)
    p <- 0
    gamma <- NULL
    vcov.method <- match.arg(vcov.method)
    if (!is.null(z)) {
        z <- as.matrix(z)
        p <- ncol(z)
        gamma <- coefficients[(len - p + 1):len]
    }
    x <- as.matrix(x)
    x1 <- cbind(1, x)
    k <- ncol(x1)
    q <- ncol(x)
    m <- (len - p)/(3 + q)
    if (round(m) != m) {
        stop("The dimension of the coefficients is incompatible with x and z. Please check the data.")
    }
    alpha <- coefficients[1:m]
    mubeta <- matrix(coefficients[(m + 1):((2 + q) * m)], nrow = q + 
        1, ncol = m)
    sigma <- coefficients[((2 + q) * m + 1):((3 + q) * m)]
    if (m == 1) {
        xz1 <- cbind(x1, z)
        I <- matrix(0, nrow = k + p + 1, ncol = k + p + 1)
        I[1:(k + p), 1:(k + p)] <- t(xz1) \%*\% xz1/sigma^2
        I[(k + p + 1), (k + p + 1)] <- n * 2/sigma^2
        if (p != 0) {
            s.1 <- c(1:k, (k + p + 1), (k + 1):(k + p))
            I <- I[s.1, s.1]
        }
        vcov <- solve(I)
        c.mat.vec <- c(rep(1, k), (1/sigma^(1/2))/2, rep(1, p))
        vcov <- diag(c.mat.vec) \%*\% vcov \%*\% diag(c.mat.vec)
    }
    else {
        sigma0 <- rep(1, m)
        mu0 <- double(m)
        an <- 1/n
        h <- 0
        epsilon <- 1e-08
        jpvt <- integer(max(k, p))
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
            y <- as.vector(y - z \%*\% gamma)
        }
        post <- matrix(out.p$post, nrow = n)
        p2 <- seq(k + 1, (k + 1) * m, by = k + 1)
        p1 <- (1:((k + 1) * m))[-p2]
        a <- diag(1/alpha[-m], nrow = m - 1, ncol = m - 1)
        a <- cbind(a, -1/alpha[m])
        abar <- a \%*\% t(post)
        xtheta <- x1 \%*\% mubeta
        Z0 <- t(t(y - xtheta)/sigma)
        f <- t(t(exp(-Z0^2/2)/sqrt(2 * pi))/sigma)
        phi <- t(t(f) * alpha)
        f0 <- rowSums(phi)
        vinv <- 1/(sigma * sigma)
        b <- t(t(Z0)/sigma)
        B <- t(vinv - t(b * b))
        c0 <- array(0, dim = c(n, m, k + 1))
        c0[, , (1:k)] <- array(tKR(x1, b), dim = c(n, m, k))
        c0[, , (k + 1)] <- -B/2
        if (vcov.method == "Hessian") {
            other.method = "OPG"
            C0 <- array(0, dim = c(n, m, k + 1, k + 1))
            x11 <- array(tKR(x1, x1), dim = c(n, k, k))
            for (i in 1:m) {
                C0[, i, (1:k), (1:k)] <- x11 * vinv[i]
                C0[, i, (1:k), k + 1] <- C0[, i, k + 1, (1:k)] <- x1 * 
                  b[, i] * vinv[i]
            }
            C0[, , k + 1, k + 1] <- t((vinv - 2 * t(B)) * vinv)/2
            Q.pi <- -abar \%*\% t(abar)
            Q.pi.theta <- matrix(0, nrow = m - 1, ncol = (k + 
                1) * m)
            for (i in 1:m) {
                zi <- a[, i] - abar
                wi <- c0[, i, ] * post[, i]
                Q.i <- colSums(tKR(wi, t(zi)))
                Q.pi.theta[, (k * (i - 1) + 1):(k * i)] <- matrix(Q.i[1:(k * 
                  (m - 1))], ncol = k)
                Q.pi.theta[, k * m + i] <- Q.i[(k * (m - 1) + 
                  1):((k + 1) * (m - 1))]
            }
            Q.theta <- matrix(0, nrow = (k + 1) * m, ncol = (k + 
                1) * m)
            for (i in 2:m) {
                for (j in 1:(i - 1)) {
                  wi <- c0[, i, ] * post[, i]
                  wj <- c0[, j, ] * post[, j]
                  Q.ij <- -colSums(tKR(wi, wj))
                  Q.theta[((k + 1) * (i - 1) + 1):((k + 1) * 
                    i), ((k + 1) * (j - 1) + 1):((k + 1) * j)] = t(matrix(Q.ij, 
                    nrow = k + 1, ncol = k + 1))
                }
            }
            Q.theta <- Q.theta + t(Q.theta)
            for (i in 1:m) {
                C.ii <- array(C0[, i, , ], dim = c(n, k + 1, 
                  k + 1))
                Q.ii.1 <- apply(C.ii * post[, i], c(2, 3), sum)
                w.ii <- tKR(c0[, i, ], c0[, i, ]) * post[, i] * 
                  (1 - post[, i])
                Q.ii.2 <- matrix(colSums(w.ii), nrow = k + 1, 
                  ncol = k + 1)
                Q.theta[((k + 1) * (i - 1) + 1):((k + 1) * i), 
                  ((k + 1) * (i - 1) + 1):((k + 1) * i)] <- -Q.ii.1 + 
                  Q.ii.2
            }
            Q.theta <- Q.theta[c(p1, p2), c(p1, p2)]
            dimI <- m - 1 + (k + 1) * m
            I <- matrix(0, nrow = dimI, ncol = dimI)
            I[1:(m - 1), 1:(m - 1)] <- -Q.pi
            I[1:(m - 1), m:dimI] <- -Q.pi.theta
            I[m:dimI, 1:(m - 1)] <- -t(Q.pi.theta)
            I[m:dimI, m:dimI] <- -Q.theta
            if (!is.null(z)) {
                dbar <- z * rowSums(post * b)
                Q.gamma.theta <- matrix(0, nrow = p, ncol = (k + 
                  1) * m)
                for (i in 1:m) {
                  C.i <- array(C0[, i, 1, ], dim = c(n, k + 1))
                  Q.i.1 <- colSums(tKR(-C.i + b[, i] * c0[, i, 
                    ], z * post[, i]))
                  Q.i.2 <- colSums(tKR(c0[, i, ] * post[, i], 
                    dbar))
                  Q.gamma.theta[, ((k + 1) * (i - 1) + 1):((k + 
                    1) * i)] <- matrix(Q.i.1 + Q.i.2, nrow = p, 
                    ncol = k + 1)
                }
                Q.gamma.theta <- Q.gamma.theta[, c(p1, p2), drop = FALSE]
                w1 <- (post * b) \%*\% t(a) - rowSums(post * b) * 
                  t(abar)
                Q.pi.gamma.0 <- colSums(tKR(w1, z))
                Q.pi.gamma <- matrix(Q.pi.gamma.0, nrow = m - 
                  1, ncol = p)
                Q.gamma <- -t(z) \%*\% (z * rowSums(post * B)) - 
                  matrix(colSums(tKR(dbar, dbar)), nrow = p, 
                    ncol = p)
                I <- cbind(I, -rbind(Q.pi.gamma, t(Q.gamma.theta)))
                I <- rbind(I, -cbind(t(Q.pi.gamma), Q.gamma.theta, 
                  Q.gamma))
            }
        }
        else {
            other.method = "Hessian"
            c0.a <- array(0, dim = c(n, m, 2))
            c0.a[, , 1] <- b
            c0.a[, , 2] <- -B/2
            score <- t(abar)
            for (j in 1:m) {
                score <- cbind(score, x1 * c0.a[, j, 1] * post[, 
                  j], c0.a[, j, 2] * post[, j])
            }
            ind <- c(1:(m - 1), p1 + m - 1, p2 + m - 1)
            score <- score[, ind]
            I <- t(score) \%*\% score
            if (!is.null(z)) {
                dbar <- z * rowSums(post * b)
                score <- cbind(score, dbar)
                I <- t(score) \%*\% score
            }
        }
        vcov <- try(solve(I))
        if (class(vcov) == "try-error" || any(diag(vcov) < 0)) {
            vcov <- matrix(NaN, nrow = (2 + k) * m - 1 + p, ncol = (2 + 
                k) * m - 1 + p)
            warning("Fisher information matrix is singular and/or the\n    variance is estimated to be negative. Consider using vcov.method=\"", 
                other.method, "\".")
        }
        c.mat.vec <- c(rep(1, m - 1 + m * k), (1/sigma^(1/2))/2, 
            rep(1, p))
        vcov <- diag(c.mat.vec) \%*\% vcov \%*\% diag(c.mat.vec)
        M.mat <- diag(len - 1)
        M.mat <- rbind(M.mat[1:(m - 1), ], c(rep(-1, m - 1), 
            rep(0, len - m)), M.mat[m:(len - 1), ])
        vcov <- M.mat \%*\% vcov \%*\% t(M.mat)
    }
    vcov
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
