#' @description Computes the critical values of the modified EM test.
#' @export
#' @title regmixCrit
#' @name regmixCrit
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gamma_1, ..., gamma_m))
#' @param z n by p matrix of regressor associated with gamma
#' @param values 3 by 1 Vector of length 3 (k = 1, 2, 3) at which the p-values are computed
#' @param parallel Determines whether package \code{doParallel} is used for calculation
#' @param cl Cluster used for parallelization; if it is \code{NULL}, the system will automatically
#' @param nrep The number of replications used to compute p-values
#' @param ninits.crit The number of initial guesses to form critical values 
#' @return A list with the following items:
#' \item{crit}{3 by 3 matrix of (0.1, 0.05, 0.01 critical values), jth row corresponding to k=j}
#' \item{pvals}{A vector of p-values at k = 1, 2, 3}
regmixCrit <- function(y, x, parlist, z = NULL, values = NULL, parallel = TRUE,
                       cl = NULL, nrep = 1000, ninits.crit = 25)
{
  # Computes the critical values of the modified EM test
  # and the estimated variance of the two-component MLE
  # Input
  #   y     : n by 1 vector of dependent variable
  #   x     : n by (q1-1) matrix of regressor not including an intercept
  #   parlist   : list including (alpha, mubeta, sigma)
  #   values (q1 by 1): values at wchich the p-values are computed
  # Output
  #   list(crit,pvals)
  #   crit = (10%, 5%, 1% critical values), pvals = p-values
  if (normalregMix.test.on) # initial values controlled by normalregMix.test.on
    set.seed(normalregMix.test.seed)

  y  <- as.vector(y)
  n  <- length(y)
  p  <- 0

  alpha   <- parlist$alpha
  mubeta  <- parlist$mubeta
  sigma   <- parlist$sigma
  gam   <- parlist$gam
  m       <- length(alpha)

  if (!is.null(z)) {
    z <- as.matrix(z)
    p <- ncol(z)
    y   <- y - z %*% gam
  }

  pvals <- NULL

  x     <- as.matrix(x)
  x1    <- cbind(1,x)
  q     <- ncol(x)

  # normalized data, n by m
  W  <- t(t(matrix(rep.int(y,m), ncol=m) - x1 %*% mubeta)/sigma)       
  f  <- t(t(exp(-W^2/2)/sqrt(2*pi))/sigma)  # pdf, n by m
  f0 <- colSums(t(f)*alpha)                 # data pdf, n by 1

  H <- hermite(W,sigma)

  w_m <- t(t(f)*alpha)/f0  # n by m matrix of w_{ji}

  if (m == 1) {
    S_alpha <- NULL
  } else {
    S_alpha <- (f[,1:(m-1)] - f[,m])/f0
  }

  S_mu    <- w_m*H[,,1]
  S_beta  <- matrix(0, nrow=n, ncol=q*m)

  for (j in 1:m) {
    S_beta[, (1+(j-1)*q):(j*q)] <- x*S_mu[, j]
  }

  S_sigma <- w_m*H[,,2]
  S_eta   <- cbind(S_alpha, S_mu, S_beta, S_sigma)

  if (!is.null(z)) {
    S_gam <- z*rowSums(S_mu)
    S_eta <- cbind(S_eta,S_gam)
  }

  n_lam <- q*(q+1)/2+2*q+2
  S_lam <- matrix(0, nrow=n, ncol=n_lam*m)
  xx    <- matrix(0, nrow=n, ncol=q*(q+1)/2)

  xx[,1:q] <- x*x

  if (q > 1) {
    t <- q + 1
    for (j in 1:(q-1)) {
      for (i in (j+1):q) {
        xx[,t] <- 2*x[,j]*x[,i]
        t <- t+1
      }
    }
  }

  for (j in 1:m) {
    w_2 <- S_sigma[,j]
    w_3 <- w_m[,j]*H[,j,3]
    w_4 <- w_m[,j]*H[,j,4]
    S_lam_1 <- cbind(w_3, x*w_2)
    S_lam_2 <- cbind(w_4, x*w_3, xx*w_2)
    S_lam[, ((j-1)*n_lam+1):(j*n_lam)] <- cbind(S_lam_1, S_lam_2)
  }

  I_eta     <- t(S_eta) %*% S_eta/n
  I_lam     <- t(S_lam) %*% S_lam/n
  I_el      <- t(S_eta) %*% S_lam/n
  if (qr(I_eta)$rank == nrow(I_eta)) {
    I_eta_lam <- I_lam - t(I_el) %*% solve(I_eta,I_el)
  } else {
    stop("The critical value cannot be computed due to singularity of some matrices.
         Please try a bootstrap version, regmixCritBoot and regmixMEMtestBoot.")
  }

  # generate u ~ N(0,I_eta_lam)
  set.seed(123456)
  e <- eigen(I_eta_lam, symmetric=TRUE)  # eigenvalue decomposition is slower than chol but more stable
  u <- t(e$vec %*% (t(e$vec) * sqrt(e$val)) %*% matrix(rnorm(nrep*n_lam*m), nrow=n_lam*m))

  q_1 <- 1+q
  q_2 <- 1+q+q*(q+1)/2

  LR <- matrix(0, nrow=nrep, ncol=m)

  if ( (parallel) && (is.null(cl)) ) {
    ncpus <- parallel::detectCores()
    cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
  }

  for (j in 1:m) {
    I_j <- I_eta_lam[((j-1)*n_lam+1):(j*n_lam),((j-1)*n_lam+1):(j*n_lam)]
    if (qr(I_j)$rank == nrow(I_j)) {
      Z_j <- u[,((j-1)*n_lam+1):(j*n_lam)] %*% solve(I_j)    # nrep by n_lam matrix
    } else {
      stop("The critical value cannot be computed due to singularity of some matrices.
           Please try a bootstrap version, regmixCritBoot and regmixMEMtestBoot.")
    }
    # Lambda_1
    V <- solve(I_j)
    n_v <- ncol(V)
    V_11 <- V[1:q_1,1:q_1]
    V_12 <- V[1:q_1,(q_1+1):n_v]
    V_21 <- t(V_12)
    V_22 <- V[(q_1+1):n_v,(q_1+1):n_v]
    Z_1 <- Z_j[,1:q_1]
    Z_2 <- Z_j[,(q_1+1):ncol(Z_j)]
    V_1_2 <- V_11 - V_12 %*% solve(V_22,V_21)
    Z_1_2 <- Z_1 - Z_2 %*% solve(V_22,V_21)
    inv_V_22 <- solve(V_22)
    Z_22 <- t(inv_V_22[1,] %*% t(Z_2))
    LR_1 <- rowSums((Z_1_2 %*% solve(V_1_2))*Z_1_2) + (1/inv_V_22[1,1])*(Z_22^2)*(Z_22<0)

    # Lambda_2
    if (parallel) {
      parallel::clusterSetRNGStream(cl, 123456)
      LR_2 <- parallel::parRapply(cl, Z_j, LR_2.comp, I_j, q, ninits.crit)
    } else {
      LR_2    <- apply(Z_j, 1, LR_2.comp, I_j, q, ninits.crit)
    }
    LR[,j]  <- apply(cbind(LR_1,LR_2), 1, max)
    }

  if (parallel) { parallel::stopCluster(cl) }

  max_EM <- apply(LR, 1, max)
  max_EM_sort <- sort(max_EM)
  qc <- floor(nrep*c(0.90,0.95,0.99))
  crit <- max_EM_sort[qc]

  if (!is.null(values)) {
    q1 <- length(values)
    pvals <- rowMeans(t(matrix(rep.int(max_EM_sort,q1),ncol=q1)) >= values)
  }

  return(list(crit = crit, pvals = pvals))

  } # end function regmixCrit


#' @description Computes the bootstrap critical values of the modified EM test.
#' @export
#' @title regmixCritBoot
#' @name regmixCritBoot
#' @param y n by 1 vector of data for y
#' @param x n by q vector of data for x
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gamma_1, ..., gamma_m))
#' @param z n by p matrix of regressor associated with gamma
#' @param values 3 by 1 Vector of length 3 (k = 1, 2, 3) at which the p-values are computed
#' @param ninits The number of initial candidates to be generated
#' @param nbtsp The number of bootstrap observations; by default, it is set to be 199
#' @param parallel Determines whether package \code{doParallel} is used for calculation
#' @param cl Cluster used for parallelization (optional)
#' @return A list with the following items:
#' \item{crit}{3 by 3 matrix of (0.1, 0.05, 0.01 critical values), jth row corresponding to k=j}
#' \item{pvals}{A vector of p-values at k = 1, 2, 3}
regmixCritBoot <- function (y, x, parlist, z = NULL, values = NULL, ninits = 100,
                            nbtsp = 199, parallel = TRUE, cl = NULL) {
  if (normalregMix.test.on) # initial values controlled by normalregMix.test.on
    set.seed(normalregMix.test.seed)

  y  <- as.vector(y)
  n  <- length(y)
  x  <- as.matrix(x)

  alpha   <- parlist$alpha
  mubeta  <- parlist$mubeta
  sigma   <- parlist$sigma
  gam   <- parlist$gam
  m       <- length(alpha)

  pvals <- NULL

  # Generate bootstrap observations
  ybset <- replicate(nbtsp, rnormregmix(n = n, x = x, alpha = alpha, mubeta = mubeta, sigma = sigma))

  if (!is.null(z)) {
    zgam <- as.matrix(z) %*% gam
    ybset <- ybset + replicate(nbtsp, as.vector(zgam))
  }

  if (parallel) {
    if (is.null(cl))
      cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    out <- foreach (j.btsp = 1:nbtsp) %dopar% {
      regmixMEMtest (ybset[,j.btsp], x = x, m = m,
                     z = z, ninits = ninits, crit.method = "none") }
    on.exit(cl)
  }
  else
    out <- apply(ybset, 2, regmixMEMtest, x = x, m = m, z = z,
                 ninits = ninits, crit.method = "none")

  emstat.b <- sapply(out, "[[", "emstat")  # 3 by nbstp matrix

  emstat.b <- t(apply(emstat.b, 1, sort))

  q <- ceiling(nbtsp*c(0.90,0.95,0.99))
  crit <- emstat.b[, q]

  if (!is.null(values)) { pvals <- rowMeans(emstat.b > values) }

  return(list(crit = crit, pvals = pvals))
}  # end function regmixCritBoot
