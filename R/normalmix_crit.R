#' @description Computes the critical values of the modified EM test.
#' @export
#' @title normalmixCrit
#' @name normalmixCrit
#' @param y n by 1 vector of data.
#' @param parlist parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam).
#' @param z n by p matrix of regressor associated with gamma.
#' @param values vector of the values of the modified EM statistic at which the p-values are computed.
#' @param nrep number of replications used to compute critical values.
#' @return A list with the following items:
#' \item{crit}{vector of critical values at the (0.1, 0.05, 0.01) level.}
#' \item{pvals}{vector of p-values corresponding to \code{values}.}
normalmixCrit <- function(y, parlist, z = NULL, values = NULL, nrep = 10000)
{
  if (normalregMixtest.env$normalregMix.test.on) # initial values controlled by normalregMix.test.on
    set.seed(normalregMixtest.env$normalregMix.test.seed)

  y <- as.vector(y)
  n <- length(y)
  p <- 0

  alpha <- parlist$alpha
  mu    <- parlist$mu
  sigma <- parlist$sigma
  gam <- parlist$gam
  m     <- length(alpha)

  if (!is.null(z)) {
    z <- as.matrix(z)
    p <- ncol(z)
    y   <- y - z %*% gam
  }

  pvals <- NULL

  if (m==1){
    crit <- qchisq(c(0.1,0.05,0.01), 2, lower.tail=F)
    if (!is.null(values))
    {
      k <- length(values)
      pvals <- pchisq(values, 2, lower.tail=F)
    }

  } else { # m>=2

    set.seed(123456)

    Z0 <- t((t(matrix(rep.int(y,m), ncol=m))-mu)/sigma)  # normalized data, n by m
    f <- t(t(exp(-Z0^2/2)/sqrt(2*pi))/sigma)      # pdf, n by m
    f0 <- colSums(t(f)*alpha)              # data pdf, n by 1

    S_alpha <- (f[, 1:(m-1)] - f[ ,m])/f0
    H <- hermite(Z0,sigma)

    S_mu    <- t(t(H[,,1]*f)*alpha)/f0    # n by m
    S_sigma <- t(t(H[,,2]*f)*alpha)/f0    # n by m
    S_gam <- rowSums(S_mu)*z

    S_lambda21 <- t(t(H[,,3]*f)*alpha)/f0  # n by m
    S_lambda22 <- t(t(H[,,4]*f)*alpha)/f0  # n by m

    S_lambda <- matrix(rbind(S_lambda21,S_lambda22),nrow=n,ncol=2*m)  # score wrt lambda (n by 2m)

    S_eta <- cbind(S_alpha,S_gam,S_mu,S_sigma)
    I11 <- t(S_eta) %*% S_eta/n      # I_eta
    I21 <- t(S_lambda) %*% S_eta/n    # I_{lambda eta}
    I22 <- t(S_lambda) %*% S_lambda/n  # I_lambda, 2*m by 2*m

    Iall <- rbind(cbind(I11,t(I21)),cbind(I21,I22))
    if (rcond(Iall) < .Machine$double.eps)
    {
      eig <- eigen(Iall, symmetric=TRUE)
      tol2 <- (1e-14)*eig$values[1]
      vals <- eig$values
      vals[vals < tol2] <- tol2
      Iall.mod <- eig$vectors %*% (vals * t(eig$vectors))
      I11 <- Iall.mod[1:(p+3*m-1),1:(p+3*m-1)]
      I21 <- Iall.mod[(p+3*m):(p+5*m-1),1:(p+3*m-1)]
      I22 <- Iall.mod[(p+3*m):(p+5*m-1),(p+3*m):(p+5*m-1)]
    }

    I221 <- I22 - I21%*%solve(I11,t(I21))  # I_{lambda.eta}, which is also (var(W_lambda.eta))^{-1}

    e <- eigen(I221, symmetric=TRUE)    # eigenvalue decomposition is slower than chol but more stable
    u <- t(e$vec %*% (t(e$vec) * sqrt(e$val)) %*% matrix(rnorm(nrep*2*m), 2*m, nrep))

    EM <- matrix(0, nrow=nrep, ncol=m)

    for (jj in 1:m) {    # computes m quadratic forms
      uu <- u[,(2*jj-1):(2*jj)]
      Ivinv <- solve(I221[(2*jj-1):(2*jj),(2*jj-1):(2*jj)])
      EM[,jj] <- rowSums((uu%*%Ivinv)*uu)
    }

    max_EM <- apply(EM,1,max)  # max of m local modified EM stats
    max_EM_sort <- sort(max_EM)
    q <- ceiling(nrep*c(0.90,0.95,0.99))
    crit <- max_EM_sort[q]

    if (!is.null(values))
    {
      k <- length(values)
      pvals <- rowMeans(t(matrix(rep.int(max_EM_sort,k),ncol=k)) >= values)
    }

  }   # end else (m>=2)
  return(list(crit = crit, pvals = pvals))
}  # end function normalmixCrit

#' @description Computes the bootstrap critical values of the modified EM test.
#' @export
#' @title normalmixCritBoot
#' @name normalmixCritBoot
#' @param y n by 1 vector of data.
#' @param parlist parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam).
#' @param z n by p matrix of regressor associated with gamma.
#' @param values vector of length 3 (k = 1, 2, 3) at which the p-values are computed.
#' @param ninits number of candidates of the initial value of the EM algorithm.
#' @param nbtsp number of bootstrap replicates. Default is 199.
#' @param parallel Determines what percentage of available cores are used, represented by a double in [0,1]. Default is 1.
#' @param cl cluster used for parallelization; if it is \code{NULL}, the system will automatically generate a cluster.
#' @param an tuning parameter used in the penalty function.
#' @return A list with the following items:
#' \item{crit}{3 by 3 matrix of critival values at the (0.1, 0.05, 0.01) level. jth row corresponding to k=j.}
#' \item{pvals}{vector of p-values at k = 1, 2, 3 corresponding to \code{values}.}
normalmixCritBoot <- function (y, parlist, z = NULL, values = NULL, ninits = 10,
                               nbtsp = 199, parallel = 1, cl = NULL, an = NULL, cn = cn) {
  if (normalregMixtest.env$normalregMix.test.on) # initial values controlled by normalregMix.test.on
    set.seed(normalregMixtest.env$normalregMix.test.seed)

  y   <- as.vector(y)
  n   <- length(y)

  alpha <- parlist$alpha
  mu    <- parlist$mu
  sigma <- parlist$sigma
  gam <- parlist$gam
  m     <- length(alpha)
  if (is.null(an)){an <- anFormula(parlist = parlist, m = m, n = n)}

  pvals <- NULL

  # Generate bootstrap observations
  ii    <- sample(m, nbtsp*n, replace = TRUE, prob = alpha)
  ybset <- rnorm(nbtsp*n, mean = mu[ii], sd = sigma[ii])
  ybset <- matrix(ybset, nrow = n, ncol = nbtsp)

  if (!is.null(z)) {
    zgam <- as.matrix(z) %*% gam
    ybset <- ybset + replicate(nbtsp, as.vector(zgam))
  }
  
  num.cores <- max(1,floor(detectCores()*parallel))
  if (num.cores > 1) {
    if (is.null(cl)) {
      cl <- makeCluster(num.cores)
      newcluster <- TRUE
    }
    out <- parCapply(cl, ybset, normalmixMEMtest, m = m, parallel = 0, 
                        z = z, an = an, cn = cn, ninits = ninits, crit.method = "none")
    if (newcluster) {
      on.exit(stopCluster(cl))
    } else {
      on.exit(cl)
    }
  }
  else
    out <- apply(ybset, 2, normalmixMEMtest, m = m, z = z, an = an, cn = cn, ninits = ninits,
                  crit.method = "none", parallel = 0)

  emstat.b <- sapply(out, "[[", "emstat")  # 3 by nbstp matrix

  emstat.b <- t(apply(emstat.b, 1, sort))

  q <- ceiling(nbtsp*c(0.90,0.95,0.99))
  crit <- emstat.b[, q]

  if (!is.null(values)) { pvals <- rowMeans(emstat.b > values) }

  return(list(crit = crit, pvals = pvals))
}  # end function normalmixCritBoot
