normalmixVcov <- function(y, coefficients, z = NULL, vcov.method = c("Hessian", "OPG"))
{
# Computes the variance-covariance matrix of the MLE of m-component normal mixture
# Input
#   y : n by 1 vector of data
#   coefficients : (alpha_1,...,alpha_m,mu_1, ...,mu_m,sigma_1, ..., sigma_m,gamma)
#   z  : n by p matrix of regressor associated with gamma
# Output
#   vcov: variance-covariance matrix

y <- as.vector(y)
n <- length(y)
len <- length(coefficients)
p <- 0
gamma  <- NULL
vcov.method <- match.arg(vcov.method)

if (!is.null(z)) {
  z <- as.matrix(z)
  p <- ncol(z)
  gamma <- coefficients[(len-p+1):len]
}

m <- (len-p)/3
if (round(m) != m) {
  stop("The dimension of the coefficients is incompatible with z. Please check the data.")
}

alpha   <- coefficients[1:m]
mu      <- coefficients[(m+1):(2*m)]
sigma   <- coefficients[(2*m+1):(3*m)]

if (m == 1) {
  if (is.null(z)) {
    I <- n*diag(c(1/sigma^2, 2/sigma^2))  # information matrix
  } else {
    z1 <- cbind(1, z)
    I <- matrix(0, nrow=p+2, ncol=p+2)
    I[1:(p+1), 1:(p+1)] <- t(z1) %*% z1/sigma^2
    I[(p+2), (p+2)] <- n*2/sigma^2
    s.1 <- c(1,(p+2), 2:(p+1))
    I <- I[s.1, s.1]
  }
  vcov  <- solve(I)
  # Because the variance is parameterized as sigma^2, we convert it to sigma
  c.mat.vec <- c(1, (1/sigma^(1/2))/2, rep(1, p))
  vcov <- diag(c.mat.vec) %*% vcov %*% diag(c.mat.vec)

} else { # end of if (m == 1)
  # m >= 2

  # Compute posterior probabilities, and adjust y if z is present
  sigma0  <- rep(1, m)  # dummy
  mu0     <- double(m)  # dummy
  an      <- 1/n  # penalty term for variance
  h       <- 0
  tau     <- 0.5
  k       <- 0
  epsilon <- 1e-08

  if (is.null(z)) {
    setting <- c(n, m, 1, 1)
    out.p <- .C("normalmixpmle", as.integer(setting), as.double(y),
                alphaset = as.double(alpha), muset = as.double(mu), sigmaset = as.double(sigma),
                as.double(sigma0), as.double(mu0), as.double(an), as.double(tau), as.integer(h),
                as.integer(k), lub = double(2*m), double(3*m), post = double(n*m),
                loglikset = double(1), penloglikset = double(1),
                notcg = integer(1), as.double(epsilon), package = "normalregMix")
  } else {
    jpvt  <- integer(p)  # pivots used in dgelsy
    setting.z <- c(n, m, p, 1, 1, jpvt)
    out.p <- .C("normalmixpmle_z", as.integer(setting.z), as.double(y), as.double(z),
                alphaset = as.double(alpha), muset = as.double(mu), sigmaset = as.double(sigma),
                gammaset = as.double(gamma),
                as.double(sigma0), as.double(mu0), as.double(an), as.double(tau), as.integer(h),
                as.integer(k), lub = double(2*m), double(3*m), post = double(n*m),
                loglikset = double(1), penloglikset = double(1),
                notcg = integer(1), as.double(epsilon), double(n*(2+p)), package = "normalregMix")
    # Adjust y
    y <- y - z %*% gamma
  }

  post <- matrix(out.p$post, nrow=n)

  p1 <- seq(1, (2*m-1), by=2) # sequence of odd numbers, 1,3,...,2*m-1
  p2 <- seq(2, (2*m), by=2)    # sequence of even numbers, 2,4,...,2*m

  # Matrices used in computing vcov

  a <- diag(1/alpha[-m], nrow=m-1, ncol=m-1)
  a <- cbind(a, -1/alpha[m])  # m-1 by m matrix of a_i's
  abar <- a %*% t(post) # m-1 by n

  Z0 <- t((t(matrix(rep.int(y, m), ncol=m))-mu)/sigma)  # normalized data, n by m
  f <- t(t(exp(-Z0^2/2)/sqrt(2*pi))/sigma)      # pdf, n by m
  phi <- t(t(f)*alpha)                # n by m
  f0 <- rowSums(phi)                  # data pdf, n by 1

  vinv <- 1/(sigma*sigma)

  b <- t(t(Z0)/sigma)  # n by m
  B <- t(vinv - t(b*b))  # n by m

  c0 <- array(0, dim=c(n, m, 2))
  c0[, , 1] <- b
  c0[, , 2] <- -B/2

  # Computes Hessian-based I
  if (vcov.method == "Hessian")  {
    other.method = "OPG"
    C0 <- array(0, dim=c(n, m, 2, 2))
    C0[, , 1, 1] <- t(matrix(vinv, nrow=m, ncol=n))  # n by m
    C0[, , 2, 1] <- C0[, , 1, 2] <- t(t(b)*vinv)    # n by m
    C0[, , 2, 2] <- t((vinv -2*t(B))*vinv)/2       # n by m

    Q.pi <- - abar %*% t(abar)  # m-1 by m-1
    Q.pi.theta <- matrix(0, nrow=m-1, ncol=2*m)  # m-1 by 2m
    for (i in 1:m){
      zi <- a[, i] - abar  # m-1 by n
      wi <- c0[, i, ]*post[, i]  # n by 2
      Q.i <- colSums(tKR(wi, t(zi)))  # 2*(m-1) vector
      # first m-1 elements correspond to mu x pi
      # second m-1 elements correspond to sigma x pi
      Q.pi.theta[, i] <- Q.i[1:(m-1)]
      Q.pi.theta[, m+i] <- Q.i[m:(2*(m-1))]
    }

    Q.theta <- matrix(0, nrow=2*m, ncol=2*m)
    for (i in 2:m){ # off-diagonal blocks
      for (j in 1:(i-1)){
        wi  <- c0[, i, ]*post[, i]
        wj  <- c0[, j, ]*post[, j]
        Q.ij <- - colSums(tKR(wi, wj))
        Q.theta[(2*i-1):(2*i), (2*j-1):(2*j)] = t(matrix(Q.ij, nrow=2, ncol=2))
      }
    }

    Q.theta <- Q.theta + t(Q.theta)
    for (i in 1:m){ # diagonal blocks
      C.ii <- array(C0[, i, , ], dim=c(n, 2, 2))
      Q.ii.1 <- apply(C.ii*post[, i], c(2, 3), sum)
      w.ii <- tKR(c0[, i, ], c0[, i, ])*post[, i]*(1-post[, i])
      Q.ii.2 <- matrix(colSums(w.ii), nrow=2, ncol=2)
      Q.theta[(2*i-1):(2*i), (2*i-1):(2*i)] <- -Q.ii.1 + Q.ii.2
    }
    # odd rows and columns of Q.theta correspond to mu
    # even rows and columns of Q.theta correspond to sigma
    Q.theta <- Q.theta[c(p1,p2), c(p1,p2)] # first block = wrt mu, second blosk = wrt sigma

    dimI <- m-1+2*m
    I <- matrix(0, nrow=dimI, ncol=dimI)
    I[1:(m-1), 1:(m-1)] <- - Q.pi
    I[1:(m-1), m:dimI]  <- - Q.pi.theta
    I[m:dimI, 1:(m-1)]  <- - t(Q.pi.theta)
    I[m:dimI, m:dimI]   <- - Q.theta

    if (!is.null(z)) {
      dbar  <-  z*rowSums(post*b) # n by p
      Q.gamma.theta <- matrix(0, nrow=p, ncol=2*m)  # p by 2*m matrix
      for (i in 1:m) {
        C.i <- array(C0[, i, 1, ], dim=c(n, 2))  # n by 2
        Q.i.1 <- colSums(tKR(-C.i+b[, i]*c0[, i, ], z*post[, i]))
        Q.i.2 <- colSums(tKR(c0[, i, ]*post[, i], dbar))  # p*2 vector
        Q.gamma.theta[, (2*i-1):(2*i)] <- matrix(Q.i.1-Q.i.2, nrow=p, ncol=2)
      }

      Q.gamma.theta <- Q.gamma.theta[, c(p1, p2),  drop=FALSE]  # p by 2*m
      w1 <- (post*b)%*%t(a) - rowSums(post*b)*t(abar)  # n by m-1
      Q.pi.gamma.0 <- colSums(tKR(w1, z))  # (m-1)*p vector
      Q.pi.gamma  <- matrix(Q.pi.gamma.0, nrow=m-1, ncol=p)
      Q.gamma     <- - t(z) %*% (z*rowSums(post*B)) -
                              matrix(colSums(tKR(dbar, dbar)), nrow=p, ncol=p)
      I <- cbind(I, -rbind(Q.pi.gamma, t(Q.gamma.theta)))
      I <- rbind(I, -cbind(t(Q.pi.gamma), Q.gamma.theta, Q.gamma))
    }  # end if (!is.null(z))

  } else  { # compute I with (vcov.method == "OPG")
    other.method = "Hessian"
    score <- t(abar)
    for (j in 1:m) { score <- cbind(score, c0[, j, ]*post[, j]) }

    ind <- c(c(1:(m-1)), p1+m-1, p2+m-1)
    score <- score[, ind]
    I <- t(score) %*% score

    if (!is.null(z)) {
      dbar  <-  z*rowSums(post*b) # n by p
      score <- cbind(score, dbar)
      I <- t(score) %*% score
    }

  } # end if (vcov.method == "OPG")

  vcov <- try(solve(I))
  if (class(vcov) == "try-error" || any(diag(vcov) <0) ) {
    vcov <- matrix(NaN, nrow = 3*m-1+p, ncol = 3*m-1+p)
    warning("Fisher information matrix is singular and/or the
    variance is estimated to be negative. Consider using vcov.method=\"",other.method,"\".")
  }

  # Because the variance is parameterized as sigma^2, we convert is to sigma

  c.mat.vec <- c(rep(1, m-1+m), (1/sigma^(1/2))/2, rep(1, p))
  vcov <- diag(c.mat.vec) %*% vcov %*% diag(c.mat.vec)

  # Add the variance of alpha_m
  len   <- length(coefficients)
  M.mat <- diag(len-1)
  M.mat <- rbind(M.mat[1:(m-1), ], c(rep(-1, m-1), rep(0, len-m)),
                 M.mat[m:(len-1), ])

  vcov     <- M.mat %*% vcov %*% t(M.mat)

}   # end else (i.e., m >= 2)

vcov

}  # end function normalmixVcov


normalmixCritBoot <- function (y, parlist, z = NULL, values = NULL, ninits = 10,
                                  nbtsp = 199, parallel = TRUE, cl = NULL) {
# Computes the bootstrap critical values of the modified EM test
# Input
#  y : n by 1 vector of data
#  coefficients : (alpha_1,...,alpha_m,mu_1, ...,mu_m,sigma_1, ..., sigma_m,gamma)
#  z  : n by p matrix of regressor associated with gamma
#  values: vector of length 3 (k = 1, 2, 3) at which the p-values are computed
# Output
#  list(crit, pvals)
#  crit = 3 by 3 matrix of (10%, 5%, 1% critical values), jth row corresponding to k=j
#  pvals = p-values at k = 1, 2, 3
SEED <- 123456
set.seed(SEED)
y   <- as.vector(y)
n   <- length(y)

alpha <- parlist$alpha
mu    <- parlist$mu
sigma <- parlist$sigma
gamma <- parlist$gamma
m     <- length(alpha)
an    <- anFormula(par = parlist, m = m, n = n)

pvals <- NULL

# Generate bootstrap observations
ii    <- sample(m, nbtsp*n, replace = TRUE, prob = alpha)
ybset <- rnorm(nbtsp*n, mean = mu[ii], sd = sigma[ii])
ybset <- matrix(ybset, nrow = n, ncol = nbtsp)

if (parallel) {
  if (is.null(cl)) {
    ncpus <- parallel::detectCores()
    cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
    parallel::clusterSetRNGStream(cl, 123456)
    out <- parallel::parCapply(cl, ybset, normalmixMEMtest, m = m, z = z, an = an,
                               ninits = ninits, crit.method = "none")
    parallel::stopCluster(cl)
  } else {
    parallel::clusterSetRNGStream(cl, 123456)
    out <- parallel::parCapply(cl, ybset, normalmixMEMtest, m = m, z = z, an = an,
                               ninits = ninits, crit.method = "none")
  }
} else {
    out <- apply(ybset, 2, normalmixMEMtest, m = m, z = z, an = an, ninits = ninits)
}

emstat.b <- sapply(out, "[[", "emstat")  # 3 by nbstp matrix

emstat.b <- t(apply(emstat.b, 1, sort))

q <- ceiling(nbtsp*c(0.90,0.95,0.99))
crit <- emstat.b[, q]

if (!is.null(values)) { pvals <- rowMeans(emstat.b > values) }

return(list(crit = crit, pvals = pvals))
}  # end function normalmixCritBoot


normalmixPMLE <- function (y, m = 2, z = NULL, vcov.method = c("Hessian", "OPG", "none"),
                           ninits = 25, epsilon = 1e-08, maxit = 2000,
                           epsilon.short = 1e-02, maxit.short = 500) {
y <- as.vector(y)
n <- length(y)
p <- 0
gamma <- NULL
ninits.short <- ninits*4*(1 + p)*m
vcov.method <- match.arg(vcov.method)

if (!is.null(z)) {
  z     <- as.matrix(z)
  p     <- ncol(z)
  if (nrow(z) != n) { stop("y and z must have the same number of rows.") }
  ls.out   <- lsfit(z, y)
  sd0    <- sqrt(mean(ls.out$residuals^2))
  jpvt  <- integer(p)  # pivots used in dgelsy
} else {
  sd0   <- sd(y) * sqrt((n - 1) / n)
}

if (m == 1) {
  if (!is.null(z)) {
    mu     <- unname(ls.out$coeff[1])
    gamma <- unname(ls.out$coeff[2:(1 + p)])
    res    <- ls.out$residuals
    sigma <- sqrt(mean(res*res))
  } else {
    mu     <- mean(y)
    sigma  <- sd0
  }

  loglik   <- - (n/2) *(1 + log(2*pi) + 2*log(sigma))
  aic      <- -2*loglik + 2*(m-1+2*m+p)
  bic      <- -2*loglik + log(n)*(m-1+2*m+p)
  penloglik <- loglik

  parlist <- list(alpha = 1, mu = mu, sigma = sigma, gamma = gamma)
  coefficients <- c(alpha = 1, mu = mu, sigma = sigma, gamma = gamma)
  postprobs <- rep(1, n)

  } else {  # m >= 2

  # generate initial values
  tmp <- normalmixPMLEinit(y = y, z = z, ninits = ninits.short, m = m)
  alphaset.s  <- tmp$alpha
  muset.s     <- tmp$mu
  sigmaset.s  <- tmp$sigma
  gammaset.s  <- tmp$gamma

  h       <- 0  # setting h=0 gives PMLE  
  k       <- 0  # k is set to 0 because this is PMLE
  tau     <- 0.5  # tau is set to 0.5 because this is PMLE
  sigma0  <- rep(sd0, m)
  mu0     <- double(m)    # dummy
  an      <- 1/n  # penalty term for variance

  # short EM
  if (is.null(z)) {
    setting <- c(n, m, ninits.short, maxit.short)  # configulation parameters
    out.short <- .C("normalmixpmle", as.integer(setting), as.double(y),
        alphaset = as.double(alphaset.s), muset = as.double(muset.s), sigmaset = as.double(sigmaset.s),
        as.double(sigma0), as.double(mu0), as.double(an), as.double(tau), as.integer(h),
        as.integer(k), lub = double(2*m), double(3*m), post = double(n*m),
        loglikset = double(ninits.short), penloglikset = double(ninits.short),
        notcg = integer(ninits.short), as.double(epsilon.short), package = "normalregMix")
  } else {
    setting.z <- c(n, m, p, ninits.short, maxit.short, jpvt)  # configulation parameters
    out.short <- .C("normalmixpmle_z", as.integer(setting.z), as.double(y), as.double(z),
        alphaset = as.double(alphaset.s), muset = as.double(muset.s), sigmaset = as.double(sigmaset.s),
        gammaset = as.double(gammaset.s),
        as.double(sigma0), as.double(mu0), as.double(an), as.double(tau), as.integer(h), 
        as.integer(k), lub = double(2*m), double(3*m), post = double(n*m),
        loglikset = double(ninits.short), penloglikset = double(ninits.short),
        notcg = integer(ninits.short), as.double(epsilon.short), double(n*(2+p)), package = "normalregMix")
  }

  penloglik.short <- out.short$penloglikset
  oo <- order(penloglik.short, decreasing = TRUE)
  oo.inits <- oo[1:ninits]

  # long EM

  alphaset  <- alphaset.s[,oo.inits]
  muset     <- muset.s[,oo.inits]
  sigmaset  <- sigmaset.s[,oo.inits]

  if (is.null(z)) {
    setting <- c(n, m, ninits, maxit)
    out <- .C("normalmixpmle", as.integer(setting), as.double(y),
        alphaset = as.double(alphaset), muset = as.double(muset), sigmaset = as.double(sigmaset),
        as.double(sigma0), as.double(mu0), as.double(an), as.double(tau), as.integer(h),
        as.integer(k), lub = double(2*m), double(3*m), post = double(n*m),
        loglikset = double(ninits), penloglikset = double(ninits),
        notcg = integer(ninits), as.double(epsilon), package = "normalregMix")
  } else {
    setting.z <- c(n, m, p, ninits, maxit, jpvt)
    gammaset  <- gammaset.s[,oo.inits]
    out <- .C("normalmixpmle_z", as.integer(setting.z), as.double(y), as.double(z),
        alphaset = as.double(alphaset), muset = as.double(muset), sigmaset = as.double(sigmaset),
        gammaset = as.double(gammaset),
        as.double(sigma0), as.double(mu0), as.double(an), as.double(tau), as.integer(h),
        as.integer(k), lub = double(2*m), double(3*m), post = double(n*m),
        loglikset = double(ninits), penloglikset = double(ninits),
        notcg = integer(ninits), as.double(epsilon), double(n*(2+p)), package = "normalregMix")
  }

  # causes an error in snow package, running on %d.
  #   if (mean(out$notcg) >= 0.9) {
  #       warning(sprintf("The EM algorithm failed to converge in %d%% of the initial values.
  #         Try increasing maxit.", 100*mean(out$notcg)))
  # }

  index     <- which.max(out$penloglikset)
  alpha     <- out$alphaset[(m*(index-1)+1):(m*index)]
  mu        <- out$muset[(m*(index-1)+1):(m*index)]
  sigma     <- out$sigmaset[(m*(index-1)+1):(m*index)]
  gamma     <- out$gammaset[(p*(index-1)+1):(p*index)]
  penloglik <- out$penloglikset[index]
  loglik    <- out$loglikset[index]

  # Compute posterior probabilities
  if (is.null(z)) {
    setting <- c(n, m, 1, 1)
    out.p <- .C("normalmixpmle", as.integer(setting), as.double(y),
        alphaset = as.double(alpha), muset = as.double(mu), sigmaset = as.double(sigma),
        as.double(sigma0), as.double(mu0), as.double(an), as.double(tau), as.integer(h),
        as.integer(k), lub = double(2*m), double(3*m), post = double(n*m),
        loglikset = double(1), penloglikset = double(1),
        notcg = integer(1), as.double(epsilon), package = "normalregMix")
  } else {
    setting.z <- c(n, m, p, 1, 1, jpvt)
    out.p <- .C("normalmixpmle_z", as.integer(setting.z), as.double(y), as.double(z),
        alphaset = as.double(alpha), muset = as.double(mu), sigmaset = as.double(sigma), gammaset = as.double(gamma),
        as.double(sigma0), as.double(mu0), as.double(an), as.double(tau), as.integer(h),
        as.integer(k), lub = double(2*m), double(3*m), post = double(n*m),
        loglikset = double(1), penloglikset = double(1),
        notcg = integer(1), as.double(epsilon), double(n*(2+p)), package = "normalregMix")
  }

  aic     <- -2*loglik + 2*(m-1+2*m+p)
  bic     <- -2*loglik + log(n)*(m-1+2*m+p)

  mu.order  <- order(mu)
  alpha     <- alpha[mu.order]
  mu        <- mu[mu.order]
  sigma     <- sigma[mu.order]

  postprobs <- matrix(out.p$post, nrow=n)
  postprobs <- postprobs[, mu.order]
  colnames(postprobs) <- c(paste("comp", ".", 1:m, sep = ""))

  parlist <- list(alpha = alpha, mu = mu, sigma = sigma, gamma = gamma)
  coefficients <- unlist(parlist)

  } # end m >= 2

if (vcov.method == "none") {
  vcov <- NULL
} else {
  vcov <- normalmixVcov(y = y, coefficients = coefficients, z = z , vcov.method = vcov.method)
}

a <- list(coefficients = coefficients, parlist = parlist, vcov = vcov, loglik = loglik,
       penloglik = penloglik, aic = aic, bic = bic, postprobs = postprobs,
       call = match.call(), m = m, label = "PMLE")

class(a) <- "normalregMix"

a

}  # end function normalmixPMLE


normalmixMEMtestSeq <- function (y, z = NULL,  maxm = 3, ninits = 10, maxit = 2000,
                          nbtsp = 199, parallel = TRUE, cl = NULL) {
# Compute the modified EM test statistic for testing H_0 of m components
# against H_1 of m+1 components for a univariate finite mixture of normals

y   <- as.vector(y)
n   <- length(y)
p   <- 0
if (!is.null(z)) {
  z <- as.matrix(z)
  p <- ncol(z)
  gamma <- matrix(0, nrow = p, ncol = maxm)
}

out   <- vector('list', length = maxm)
aic    <- bic <- double(maxm)
pvals   <- emstat <- matrix(0, nrow = maxm-1, ncol = 3)
loglik  <- penloglik <- double(maxm)

alpha   <- mu <- sigma <- matrix(0, nrow = maxm, ncol = maxm)

# Test H_0:m=1, H_0:m=2, ...

for (m in 1:maxm){
  par0   <- normalmixPMLE(y = y, m = m, z = z, vcov.method = "none",
                            ninits = ninits, maxit = maxit)
  loglik[m] <- loglik0 <- par0$loglik
  penloglik[m] <- penloglik0 <- par0$penloglik
  aic[m]  <- par0$aic
  bic[m]  <- par0$bic

  parlist <- par0$parlist
  alpha0  <- parlist$alpha
  mu0     <- parlist$mu
  sigma0  <- parlist$sigma
  gamma0  <- parlist$gamma

  alpha[,m] <- c(alpha0, double(maxm - m))
  mu[,m]     <- c(mu0, double(maxm - m))
  sigma[,m] <- c(sigma0, double(maxm - m))

  cat(sprintf("%d-component model estimate:\n",m))
  tab = as.matrix(rbind(alpha0, mu0, sigma0))
  rownames(tab) <- c("alpha", "mu", "sigma")
  colnames(tab) <- c(paste("comp", ".", 1:m, sep = ""))
  print(tab, digits = 4)

  if (!is.null(z)){
    gamma[, m] <- gamma0
    cat("gamma =", gamma0,"\n")
  }
  cat(sprintf("\nAIC, BIC, and loglike of 1 to %.i", m), "component models \n")
  cat(c("AIC    =", sprintf(' %.2f', aic[1:m])), "\n")
  cat(c("BIC    =", sprintf(' %.2f', bic[1:m])), "\n")
  cat(c("loglik =", sprintf('%.2f', loglik[1:m])), "\n\n")

  if (m < maxm){

    cat(sprintf("Testing the null hypothesis of %d components\n", m))
  
#     print(par0$parlist)
    an    <- anFormula(par = parlist, m = m, n = n)
#     print(an)
    par1  <- normalmixMaxPhi(y = y, par0 = parlist, z = z, an = an, ninits = ninits, maxit = maxit)
#     print(loglik0)
#     print(par1)
    emstat.m  <- 2*(par1$loglik - loglik0)
  
    cat(c("modified EM-test statitic ", sprintf('%.3f',emstat.m)),"\n")
    if (m <=3 ) {
      em.out <- normalmixCrit(y=y, parlist=parlist, z=z, values = emstat.m)
      cat(c("asymptotic p-values       ", sprintf('%.3f',em.out$pvals)),"\n \n")
    } else {
      em.out <- normalmixCritBoot(y=y, parlist=parlist, z=z, values = emstat.m,
                    ninits = ninits, nbtsp = nbtsp, parallel = parallel, cl = cl)
      cat(c("bootstrap p-values        ", sprintf('%.3f',em.out$pvals)),"\n \n")
    }
    # noncg.rate[m]   <- par1$noncg.rate
    pvals[m,]     <- em.out$pvals
    emstat[m,]    <- emstat.m
  }
}
a = list(alpha = alpha, mu = mu, sigma = sigma, gamma = gamma, emstat = emstat, pvals = pvals, aic = aic, bic = bic, loglik = loglik, penloglik = penloglik)

}  # end normalmixMEMtestSeq


normalmixMEMtest <- function (y, m = 2, z = NULL, an = NULL, tauset = c(0.1,0.3,0.5), 
                              ninits = 10,
                              crit.method = c("asy", "boot", "none"), nbtsp = 199,
                              cl = NULL, 
                              parallel.method = c("none", "do", "snow")) {
# Compute the modified EM test statistic for testing H_0 of m components
# against H_1 of m+1 components for a univariate finite mixture of normals

y     <- as.vector(y)
n     <- length(y)
if (!is.null(z)) {z <- as.matrix(z)}
crit.method <- match.arg(crit.method)
parallel.method <- match.arg(parallel.method)

par0    <- normalmixPMLE(y=y, m=m, z=z, vcov.method = "none", ninits=ninits)
loglik0 <- par0$loglik

if (is.null(an)){ an <- anFormula(par = par0$parlist, m = m, n = n) }

par1  <- normalmixMaxPhi(y=y, par0=par0$parlist, z=z, 
                         an=an, tauset = tauset, ninits=ninits,
                         parallel.method = parallel.method, cl = cl)


emstat  <- 2*(par1$penloglik - loglik0)
# emstat  <- 2*(par1$loglik - loglik0)
# print(par1)

if (crit.method == "asy"){
  result  <- normalmixCrit(y=y, parlist=par0$parlist, z=z, values=emstat)
} else if (crit.method == "boot") {
  result  <- normalmixCritBoot(y=y, parlist= par0$parlist, z=z, values=emstat,
                            ninits=ninits, nbtsp=nbtsp, (parallel.method != "none"), cl=cl)
} else {
  result <- list()
  result$crit <- result$pvals <- rep(NA,3)
}

a <- list(emstat = emstat, pvals = result$pvals, crit = result$crit, crit.method = crit.method,
          parlist = par0$parlist, call = match.call(), m = m, label = "MEMtest")

class(a) <- "normalregMix"

a

}  # end normalmixMEMtest



normalmixMaxPhi <- function (y, par0, z = NULL, an, tauset = c(0.1,0.3,0.5), 
                             ninits = 10, epsilon.short = 1e-02, epsilon = 1e-08,
                             maxit.short = 500, maxit = 2000, 
                             verb = FALSE, 
                             parallel.method = c("none", "do", "snow"),
                             cl = NULL) {
# Given a parameter estimate of an m component model and tuning paramter an,
# maximize the objective function for computing the modified EM test statistic
# for testing H_0 of m components against H_1 of m+1 for a univariate normal finite mixture

warn  <- options(warn=-1) # Turn off warnings
parallel.method <- match.arg(parallel.method)

p <- 0
m <- length(par0$alpha)

# jpvt must be specified to run SNOW package
jpvt <- NULL
if (!is.null(z)) {
  p <- ncol(z)
  jpvt   <- integer(p)    # pivots used in dgelsy
}

ninits.short <- ninits*4*(1+p)*m



loglik.all <- matrix(0,nrow=m*length(tauset),ncol=3)
penloglik.all <- matrix(0,nrow=m*length(tauset),ncol=3)
# Remember, you can always take phi1 back,
# phi1 <- vector('list',m*length(tauset))

if (parallel.method == "do") {
  # TODO: Assign or check with cl, if possible.
  registerDoParallel()
  results <- foreach (t = 1:length(tauset),
                      .export = 'normalmixMaxPhiStep', .combine = c)  %:%
    foreach (h = 1:m) %dopar% {
      normalmixMaxPhiStep (c(h, tauset[t]), y, par0, z = NULL, p, jpvt,
                           an,
                           ninits, ninits.short,
                           epsilon.short, epsilon,
                           maxit.short, maxit,
                           verb) }
  loglik.all <- t(sapply(results, "[[", "loglik"))
  penloglik.all <- t(sapply(results, "[[", "penloglik"))
}
else if (parallel.method == "snow") {
  if (is.null(cl))
    cl <- snow::makeCluster(detectCores(), type = "SOCK")
  
  # causes an error in snow package, running on %d.
  htaupairs <- expand.grid(h = seq(1,m), tau = tauset)
  results <- snow::parApply(cl, htaupairs, 1, normalmixMaxPhiStep,
                        y, par0, z, p, jpvt,
                        an,
                        ninits, ninits.short,
                        epsilon.short, epsilon,
                        maxit.short, maxit,
                        verb)
  snow::stopCluster(cl)
  loglik.all <- t(sapply(results, "[[", "loglik"))
  penloglik.all <- t(sapply(results, "[[", "penloglik")) 
}
else
  for (h in 1:m) 
    for (t in 1:length(tauset)) {
      rowindex <- (t-1)*m + h
      tau <- tauset[t]
      result <- normalmixMaxPhiStep(c(h, tau), y, par0, z = NULL, p, jpvt, 
                                    an,
                                    ninits, ninits.short,
                                    epsilon.short, epsilon,
                                    maxit.short, maxit,
                                    verb)
      # phi1[[rowindex]] <- result$phi1
      loglik.all[rowindex,] <- result$loglik
      penloglik.all[rowindex,] <- result$penloglik
    }
loglik <- apply(loglik.all, 2, max)  # 3 by 1 vector
penloglik <- apply(penloglik.all, 2, max)  # 3 by 1 vector

out <- list(loglik = loglik, penloglik = penloglik)

out

}  # end normalmixMaxPhi

normalmixMaxPhiStep <- function (htaupair, y, par0, z = NULL, p, jpvt, 
                                 an,
                                 ninits, ninits.short,
                                 epsilon.short, epsilon,
                                 maxit.short, maxit,
                                 verb)
{
  # TODO: Beautify and move some parameters to the method
  # that calls this.
  alpha0 <- par0$alpha
  
  m      <- length(alpha0)
  m1     <- m+1
  k      <- 1
  n      <- length(y)
  h      <- as.numeric(htaupair[1])
  tau    <- as.numeric(htaupair[2])
  
  mu0    <- par0$mu
  mu0h   <- c(0,mu0,0)        # m+1 by 1
  sigma0 <- par0$sigma
  sigma0h<- c(sigma0[1:h],sigma0[h:m])
  gamma0 <- par0$gamma

  # generate initial values
  tmp <- normalmixPhiInit(y = y, par = par0, z = z, h=h, tau = tau, ninits = ninits.short)
  # tmp$alpha, tmp$sigma: m+1 by ninits matrix, tmp$mu: p by m+1 by ninits array
  # tmp$gamma: p by ninits matrix
  
  alphaset.s   <- tmp$alpha
  muset.s   <- tmp$mu
  sigmaset.s   <- tmp$sigma
  gammaset.s  <- tmp$gamma

  
  if (is.null(z)) {
    setting <- c(n,m1,ninits.short,maxit.short)  # configulation parameters
    out.short <- .C("normalmixpmle", as.integer(setting), as.double(y),
                    alphaset = as.double(alphaset.s), muset = as.double(muset.s), sigmaset = as.double(sigmaset.s),
                    as.double(sigma0h), as.double(mu0h), as.double(an), as.double(tau), as.integer(h),
                    as.integer(k), lub = double(2*m1), double(3*m1), post = double(n*m1),
                    loglikset = double(ninits.short), penloglikset = double(ninits.short),
                    notcg = integer(ninits.short), as.double(epsilon.short), package = "normalregMix")
  } else {
    setting.z <- c(n,m1,p,ninits.short,maxit.short,jpvt)  # configulation parameters
    out.short <- .C("normalmixpmle_z", as.integer(setting.z), as.double(y), as.double(z),
                    alphaset = as.double(alphaset.s), muset = as.double(muset.s), sigmaset = as.double(sigmaset.s), gammaset = as.double(gammaset.s),
                    as.double(sigma0h), as.double(mu0h), as.double(an), as.double(tau), as.integer(h),
                    as.integer(k), lub = double(2*m1), double(3*m1), post = double(n*m1),
                    loglikset = double(ninits.short), penloglikset = double(ninits.short),
                    notcg = integer(ninits.short), as.double(epsilon.short), double(n*(2+p)), package = "normalregMix")
  }
  
  if (verb && any(out.short$notcg)) {
    cat(sprintf("non-convergence rate at short-EM = %.3f\n",mean(out.short$notcg)))
  }
  
  penloglik.short <- out.short$penloglikset    # extract the 4th argument = penloglik
  oo <- order(penloglik.short, decreasing = TRUE)
  oo.inits <- oo[1:ninits]
  
  # long EM
  
  alphaset   <- alphaset.s[,oo.inits]
  muset      <- muset.s[,oo.inits]
  sigmaset   <- sigmaset.s[,oo.inits]
  
  if (is.null(z)) {
    setting <- c(n,m1,ninits,maxit)
    out <- .C("normalmixpmle", as.integer(setting), as.double(y),
              alphaset = as.double(alphaset), muset = as.double(muset), sigmaset = as.double(sigmaset),
              as.double(sigma0h), as.double(mu0h), as.double(an), as.double(tau), as.integer(h),
              as.integer(k), lub = double(2*m1), double(3*m1), post = double(n*m1),
              loglikset = double(ninits), penloglikset = double(ninits),
              notcg = integer(ninits), as.double(epsilon), package = "normalregMix")
  } else {
    gammaset  <- gammaset.s[,oo.inits]
    setting.z <- c(n,m1,p,ninits,maxit,jpvt)
    out <- .C("normalmixpmle_z", as.integer(setting.z), as.double(y), as.double(z),
              alphaset = as.double(alphaset), muset = as.double(muset), sigmaset = as.double(sigmaset), gammaset = as.double(gammaset),
              as.double(sigma0h), as.double(mu0h), as.double(an), as.double(tau), as.integer(h),
              as.integer(k), lub = double(2*m1), double(3*m1), post = double(n*m1),
              loglikset = double(ninits), penloglikset = double(ninits),
              notcg = integer(ninits), as.double(epsilon), double(n*(2+p)), package = "normalregMix")
  }
  # causes an error in snow package, running on %d.
  # if (mean(out$notcg) >= 0.9) {
  #   warning(sprintf("The EM algorithm failed to converge in %d%% of the
  #                   initial values at h = %d. Please increase maxit.", 100*mean(out$notcg), h))
  # }
  
  index       <- which.max(out$penloglikset)
  alpha       <- out$alphaset[(m1*(index-1)+1):(m1*index)]
  mu          <- out$muset[(m1*(index-1)+1):(m1*index)]
  sigma       <- out$sigmaset[(m1*(index-1)+1):(m1*index)]
  gamma       <- out$gammaset[(p*(index-1)+1):(p*index)]
  penloglik   <- out$penloglikset[index]
  loglik      <- out$loglikset[index]
  noncg.rate  <- mean(out$notcg)
  
  #    print(c(h,k,tau,alpha,mu,sigma))
  #   print(loglik)
  #   print(penloglik)
  
  mu.order  <- order(mu)
  alpha     <- alpha[mu.order]
  mu        <- mu[mu.order]
  sigma     <- sigma[mu.order]
  
  a0 = list(alpha = alpha, mu = mu, sigma = sigma, gamma = gamma,
            penloglik = penloglik, loglik=loglik)
  
  # update parameters
  a     <- normalmixPhi2(y, z, a0, sigma0 = sigma0, h=h, tau=tau, an=an)
  a[[1]]  <- list(alpha = alpha, mu = mu, sigma = sigma, gamma = gamma,
                  penloglik = a0$penloglik, loglik=a0$loglik)
  
  phi1 <- c(a[[1]], h = h)
  loglik <- sapply(a, "[[", "loglik")    # extract loglik at k=1,2,3
  penloglik <- sapply(a, "[[", "penloglik")    # extract penloglik at k=1,2,3
  
  return (list(phi1 = phi1, loglik = loglik, penloglik = penloglik))
}

normalmixPhi2 <- function (y, z = NULL, par0, sigma0, h, tau, an) {
# Starting from a parameter estimate, take two EM steps to compute the
# local modified EM test statistic for testing H_0 of m components
# against H_1 of m+1 at K=2,3 for a univariate normal finite mixture
# In input,
#  y    : n by 1 vector of dependent variable
#  z    : n by p matrix of regressor associated with gamma

n     <- length(y)
alpha <- par0$alpha
mu    <- par0$mu
sigma <- par0$sigma
gamma <- par0$gamma
m1    <- length(alpha) # m1 = m+1
m     <- m1-1
sigma0h <- c(sigma0[1:h],sigma0[h:m])
a     <- vector('list',3)
mu0   <- double(m1)
h0    <- 0

p      <- 0
if (!is.null(z)) {
  p <- ncol(z)
  jpvt <- integer(p)
  setting.z <- c(n,m1,p,1,1,jpvt)
}

# setting <- c(n,m1,1,2)
setting <- c(n,m1,1,1)

    for (k in 2:3) {
    a[[k]] <- list(penloglik = -Inf, loglik = -Inf)
    }

    for (k in 2:3) {
    # Two EM steps
  if (is.null(z)) {
    out <- .C("normalmixpmle", as.integer(setting), as.double(y),
              alpha = as.double(alpha), mu = as.double(mu), sigma = as.double(sigma),
              as.double(sigma0h), as.double(mu0), as.double(an), as.double(tau), as.integer(h),
              as.integer(k), lub = double(2*m1), double(3*m1), post = double(n*m1),
              loglik = double(1), penloglik = double(1),
              notcg = integer(1), tol = as.double(1e-8), package = "normalregMix")
  } else {
    out <- .C("normalmixpmle_z", as.integer(setting.z), as.double(y), as.double(z),
              alpha = as.double(alpha), mu = as.double(mu), sigma = as.double(sigma), gamma = as.double(gamma),
              as.double(sigma0h), as.double(mu0), as.double(an), as.double(tau), as.integer(h),
              as.integer(k), lub = double(2*m1), double(3*m1), post = double(n*m1),
              loglik = double(1), penloglik = double(1),
              notcg = integer(1), tol = as.double(1e-8), double(n*(2+p)), package = "normalregMix")
  }

       # Check singularity: if singular, break from the loop
    if ( any(out$sigma < 1e-06) || any(out$alpha < 1e-06) || is.na(sum(out$alpha)) ) {
      break
    }

    alpha    <- out$alpha
    mu      <- out$mu
    sigma    <- out$sigma
    gamma    <- out$gamma
    loglik    <- out$loglik
    penloglik   <- out$penloglik

    a[[k]] <- list(alpha = alpha, mu = mu, sigma = sigma, gamma = gamma, penloglik = penloglik, loglik = loglik)

      }

  a

}  # end function normalmixPhi2


normalmixPhiInit <- function (y, par, z = NULL, h, tau, ninits = 1)
{
# Generate initial values used by the modified EM test for mixture of normals
# input is
#  y       : data
#  par     : list of alpha, mu, sigma from an m-component model
#  z       : p-dimensional covariate whose coefficient is common across components
#  h       : h in the EM test
#  tau     : parameter used to split the h-th component
#  ninits  : number of initial values to be generated
# output is a list that contains
#  mu     : m+1 by ninits matrix
#  sigma  : m+1 by ninits matrix
#  alpha  : m+1 by ninits matrix
#  gamma  : p by ninits matrix

n     <- length(y)
p     <- ncol(z)
gamma <- NULL

mu0      <- par$mu
sigma0   <- par$sigma
alpha0   <- par$alpha
gamma0  <- par$gamma
m       <- length(alpha0)

if (!is.null(z)){
  y    <- y - z %*% gamma0
  gamma <- matrix(runif(p*ninits,min=0.5,max=1.5),nrow=p)*gamma0
}

if (m>=2){
  mid <- (mu0[1:(m-1)]+mu0[2:m])/2  # m-1 by 1
  lb0 <- c(min(y),mid)        # m by 1
  lb  <- c(lb0[1:h],lb0[h:m])      # m+1 by 1
  ub0 <- c(mid,max(y))        # m by 1
  ub  <- c(ub0[1:h],ub0[h:m])      # m+1 by 1
} else {
  lb  <- c(min(y),min(y))
  ub  <- c(max(y),max(y))
}

mu <- matrix(runif((m+1)*ninits,min=lb,max=ub),nrow=m+1)

sigma.hyp <- c(sigma0[1:h],sigma0[h:m])  # m+1 by 1
sigma <- matrix(runif((m+1)*ninits,min=sigma.hyp*0.25,max=sigma.hyp*2),nrow=m+1)

alpha.hyp <- c(alpha0[1:h],alpha0[h:m])  # m+1 by 1
alpha.hyp[h:(h+1)] <- c(alpha.hyp[h]*tau,alpha.hyp[h+1]*(1-tau))
alpha <- matrix(rep.int(alpha.hyp,ninits),nrow=m+1)

list(alpha = alpha, mu = mu, sigma = sigma, gamma = gamma)

}  # end function normalmixPhiInit


normalmixCrit <- function(y, parlist, z = NULL, values = NULL, nrep = 10000)
{
# Computes the critical values of the modified EM test
# Input
#  y (n by 1): data
#  z : covariates whose coefficients are common across components
#  coefficients:
#  values (k by 1): values at wchich the p-values are computed
# Output
#  list(crit, pvals)
#  crit = (10%, 5%, 1% critical values), pvals = p-values
y <- as.vector(y)
n <- length(y)
p <- 0

alpha <- parlist$alpha
mu    <- parlist$mu
sigma <- parlist$sigma
gamma <- parlist$gamma
m     <- length(alpha)

if (!is.null(z)) {
  z <- as.matrix(z)
  p <- ncol(z)
  y   <- y - z %*% gamma
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
S_gamma <- rowSums(S_mu)*z

S_lambda21 <- t(t(H[,,3]*f)*alpha)/f0  # n by m
S_lambda22 <- t(t(H[,,4]*f)*alpha)/f0  # n by m

S_lambda <- matrix(rbind(S_lambda21,S_lambda22),nrow=n,ncol=2*m)  # score wrt lambda (n by 2m)

S_eta <- cbind(S_alpha,S_gamma,S_mu,S_sigma)
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


normalmixPMLEinit <- function (y, z = NULL, ninits = 1, m = 2)
{
# Generate initial values used by the PMLE of mixture of normals
# input
#  y    : data
#  z    : p-dimensional covariate whose coefficient is common across components
#  ninits  : number of initial values to be generated
#  m    : number of components
# output is a list that contains
#  alpha  : m by ninits matrix
#  mu     : m by ninits matrix
#  sigma  : m by ninits matrix
#  gamma  : p by ninits matrix

n <- length(y)
p <- ncol(z)
gamma <- NULL
if (!is.null(z)){
  out     <- lsfit(z,y)
  gamma0  <- out$coef[-1]
  gamma   <- matrix(runif(p*ninits, min=0.5, max=1.5), nrow=p) * gamma0
  y       <- out$residuals + out$coef[1]
}

alpha <- matrix(runif(m * ninits), nrow=m)
alpha <- t(t(alpha) / colSums(alpha))
mu    <- matrix(runif(m * ninits, min = min(y), max = max(y)), nrow=m)
sigma <- matrix(runif(m * ninits, min = 0.01, max = 2)*sd(y), nrow=m)

list(alpha = alpha, mu = mu, sigma = sigma, gamma = gamma)

}  # end function normalmixPMLEinit

rnormregmix <- function (n, x = NULL, alpha, mubeta, sigma) {
# Generates mixed normal random variables with regressor x
# Input
#  n : number of observations
#   x : (n by k-1) matrix NOT including a constant
#   alpha  : m-vector
#  mubeta  : k by m matrix
#  sigma  : m-vector
# Output
#  y : n by 1 vector

m     <- length(alpha)
mubeta   <- matrix(mubeta, ncol=m)

if (!is.null(x)){
    x <- as.matrix(x)
  if (nrow(x) != n) { stop("y and x must have the same number of rows.") }
  x1   <- cbind(1,x)
  ii   <- sample(m, n, replace=TRUE, prob=alpha)
  y   <- rnorm(n, mean = rowSums(x1*t(mubeta[, ii])), sd = sigma[ii])
} else {
  ii   <- sample(m, n, replace=TRUE, prob=alpha)
  y   <- rnorm(n, mean = mubeta[, ii], sd = sigma[ii])
}

y

}  # end function rnormregmix

omega.ji <- function(alpi, mui, sigi, alpj, muj, sigj)
# Computes omega_{j|i} defined in (2.1) of Maitra and Melnykov
{

if(sigi == sigj)
{
delta <- abs(mui - muj) / sigi
out <- pnorm(-delta / 2 + log(alpj / alpi) / delta, 0, 1)
}
else
{
ncp <- (mui - muj)*sigi / (sigi^2 - sigj^2)

value <- sigj^2 * (mui - muj)^2 / (sigj^2 - sigi^2)^2
      - sigj^2 / (sigi^2 - sigj^2) * log(alpi^2 * sigj^2 / alpj^2 / sigi^2 )
sqrt.value <- sqrt(max(value,0))

ind <- as.numeric(sigi<sigj)
out <- ind + (-1)^ind * (pnorm(sqrt.value - ncp,0,1) - pnorm(-sqrt.value - ncp,0,1))
}

return(out)

}  # end function omega.ji


omega.12 <- function(par)
# Computes omega_{12} for testing H_0:m=2 against H_1:m=3
{
alp1 <- par$alpha[1]
alp2 <- par$alpha[2]

mu1 <- par$mu[1]
mu2 <- par$mu[2]

sig1 <- par$sigma[1]
sig2 <- par$sigma[2]

part1 <- omega.ji(alp1, mu1, sig1, alp2, mu2, sig2)
part2 <- omega.ji(alp2, mu2, sig2, alp1, mu1, sig1)

return((part1 + part2) / 2)
}  # end function omega.12


omega.123 <- function(par)
{
alp1 <- par$alpha[1]
alp2 <- par$alpha[2]
alp3 <- par$alpha[3]

mu1 <- par$mu[1]
mu2 <- par$mu[2]
mu3 <- par$mu[3]

sig1 <- par$sigma[1]
sig2 <- par$sigma[2]
sig3 <- par$sigma[3]

part1 <- omega.ji(alp1, mu1, sig1, alp2, mu2, sig2)
part2 <- omega.ji(alp2, mu2, sig2, alp1, mu1, sig1)
w12 <- (part1 + part2)/2

part3 <- omega.ji(alp2, mu2, sig2, alp3, mu3, sig3)
part4 <- omega.ji(alp3, mu3, sig3, alp2, mu2, sig2)
w23 <- (part3 + part4)/2

return(c(w12, w23))

}  # end function omega.123


anFormula <- function(par, m, n)
# Computes a_n for testing H_0 of m components
# against H_1 of m+1 components
{
if (m == 1) {
#   an <- 1.0
  an <- 0.40
#   an <- 0.25
}
else if (m == 2) {
  omega <- omega.12(par)
  omega <- pmin(pmax(omega, 1e-16), 1-1e-16)  # an becomes NaN if omega[j]=0 or 1
  x <- exp(-1.747 - 0.297 * log(omega / (1 - omega)) - 98.35/n)  # maxa=1
  an <- 0.9 * x / (1 + x)
#   x <- exp(-1.642 - 0.434 * log(omega / (1 - omega)) - 101.80/n)  # maxa=2
#   an <- 1.8 * x / (1 + x)
}
else if (m == 3) {
  omega <- omega.123(par)
  omega <- pmin(pmax(omega, 1e-16), 1-1e-16)  # an becomes NaN if omega[j]=0 or 1
  t_omega <- (omega[1] * omega[2])/(1 - omega[1])/(1 - omega[2])
  x <- exp(-1.731 - 0.162 * log(t_omega) - 144.64/n)
  an <- 0.80 * x / (1 + x)
#   x <- exp(-1.678 - 0.232 * log(t_omega) - 175.50/n)
#   an <- 1.5 * x / (1 + x)
}
else {
  an <- 1.0
}

}  # end function anFormula


coef.to.list <- function(coefficients, z = NULL) {
# 　Convert coefficients to list
len     <- length(coefficients)
p       <- 0
gamma   <- NULL

if (!is.null(z)) {
  z <- as.matrix(z)
  p <- ncol(z)
  gamma <- coefficients[(len-p+1):len]
}

m <- (len-p)/3
if (round(m) != m) {
  stop("The dimension of the coefficients is incompatible with z. Please check the data.")
}

param   <- matrix(coefficients[1:(len-p)], nrow=m, ncol=3)
alpha   <- param[, 1]
mu      <- param[, 2]
sigma   <- param[, 3]

a = list(alpha = alpha, mu = mu, sigma = sigma, gamma = gamma)

a

}
