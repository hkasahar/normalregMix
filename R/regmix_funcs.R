#' @description Compute ordinary & penalized log-likelihood ratio resulting from 
#' MEM algorithm at k=1,2,3.
#' @export
#' @title regmixMaxPhi
#' @name regmixMaxPhi
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gam
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gam_1, ..., gam_m))
#' @param z n by p matrix of regressor associated with gam
#' @param an a term used for penalty function
#' @param tauset A set of initial tau value candidates
#' @param ninits The number of randomly drawn initial values.
#' @param epsilon.short The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param epsilon The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}.
#' @param maxit.short The maximum number of iterations in short EM.
#' @param maxit The maximum number of iterations.
#' @param verb Determines whether to print a message if an error occurs.
#' @return A list with items:
#' \item{loglik}{Log-likelihood resulting from MEM algorithm at k=1,2,3.}
#' \item{penloglik}{Penalized log-likelihood resulting from MEM algorithm at k=1,2,3.}
regmixMaxPhi <- function (y, x, parlist, z = NULL, an, tauset = c(0.1,0.3,0.5),
                          ninits = 10,
                          epsilon.short = 1e-02, epsilon = 1e-08,
                          maxit.short = 500, maxit = 2000,
                          verb = FALSE) {
  # Given a parameter estiamte of an m component model and tuning paramter an,
  # maximize the objective function for computing the modified EM test statistic
  # for testing H_0 of m components against H_1 of m+1 for a univariate finite mixture of normals

  warn  <- options(warn=-1) # Turn off warnings

  q1 <- ncol(x) + 1
  m <- length(parlist$alpha)

  if (!is.null(z)) {
    z     <- as.matrix(z)
    p     <- ncol(z)
  } 
  else
    p <- 0

  ninits.short <- ninits*10*(q1+p)*m

  htau <- t(as.matrix(expand.grid(c(1:m), tauset)))
  results <- apply(htau,2,regmixPhiStep, y, x, parlist, z, p, an,
                    ninits, ninits.short, epsilon.short, epsilon,
                    maxit.short, maxit, verb)
  loglik.all <- t(sapply(results, "[[", "loglik"))
  penloglik.all <- t(sapply(results, "[[", "penloglik"))
  coefficient.all <- t(sapply(results, "[[", "coefficient"))
    
  loglik <- apply(loglik.all, 2, max)  # 3 by 1 vector
  penloglik <- apply(penloglik.all, 2, max)  # 3 by 1 vector
  index <- which.max(loglik.all[ ,3]) # a par (h,m) that gives the highest likelihood at k=3
  coefficient <- as.vector(coefficient.all[index,])

  out <- list(coefficient = coefficient, loglik = loglik, penloglik = penloglik)

  out

}  # end regmixMaxPhi

#' @description Given a pair of h and tau and data, compute ordinary &
#' penalized log-likelihood ratio resulting from MEM algorithm at k=1,2,3, 
#' tailored for parallelization.
#' @export
#' @title regmixPhiStep
#' @name regmixPhiStep
#' @param htaupair A set of h and tau
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gam
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gam_1, ..., gam_m))
#' @param z n by p matrix of regressor associated with gam
#' @param p Dimension of z
#' @param an a term used for penalty function
#' @param ninits The number of randomly drawn initial values.
#' @param ninits.short The number of candidates used to generate an initial phi, in short MEM
#' @param epsilon.short The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param epsilon The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}.
#' @param maxit.short The maximum number of iterations in short EM.
#' @param maxit The maximum number of iterations.
#' @param verb Determines whether to print a message if an error occurs.
#' @return A list of coefficients, log-likelihood, and penalized log-likelihood resulting from MEM algorithm.
regmixPhiStep <- function (htaupair, y, x, parlist, z = NULL, p,
                          an,
                          ninits, ninits.short,
                          epsilon.short, epsilon,
                          maxit.short, maxit,
                          verb)
{
  alpha0 <- parlist$alpha

  x1    <- cbind(1,x)
  q1     <- ncol(x1)
  m      <- length(alpha0)
  m1     <- m+1
  n      <- length(y)
  h      <- as.numeric(htaupair[1])
  tau    <- as.numeric(htaupair[2])
  k <- 1

  mubeta0 <- parlist$mubeta
  mu0h <- c(-1e+10,mubeta0[1,],1e+10)        # m+1 by 1
  sigma0 <- parlist$sigma
  sigma0h <- c(sigma0[1:h],sigma0[h:m])  # m+1 by 1
  gam0 <- parlist$gam
  
  if (is.null(z)) {
    ztilde <- matrix(0) # dummy
    gam <- NULL
  } else
    ztilde <- as.matrix(z)

  # generate initial values
  tmp <- regmixPhiInit(y = y, x = x, z = z, parlist=parlist, h=h, tau, ninits = ninits.short)

  # short EM
  b0 <- as.matrix(rbind( tmp$alpha, tmp$mubeta, tmp$sigma, tmp$gam ))
  out.short <- cppRegmixPMLE(b0, y, x, ztilde, mu0h, sigma0h, m1, p, an, maxit.short,
                             ninits.short, epsilon.short, tau, h, k)
  # long EM
  components <- order(out.short$penloglikset, decreasing = TRUE)[1:ninits]
  b1 <- as.matrix(b0[ ,components]) # b0 has been updated
  out <- cppRegmixPMLE(b1, y, x, ztilde, mu0h, sigma0h, m1, p, an, maxit, ninits, epsilon, tau, h, k)

  index     <- which.max(out$penloglikset)
  alpha <- b1[1:m1,index] # b0 has been updated
  mubeta <- matrix(b1[(1+m1):((q1+1)*m1),index],nrow=q1,ncol=m1)
  sigma <- b1[(1+(q1+1)*m1):((q1+2)*m1),index]
  if (!is.null(z)) {
    gam     <- b1[((q1+2)*m1+1):((q1+2)*m1+p),index]
  }
  mu.order  <- order(mubeta[1,])
  alpha     <- alpha[mu.order]
  mubeta    <- mubeta[ ,mu.order]
  sigma     <- sigma[mu.order]
  sigma0h <- sigma0h[mu.order]
  b <- as.matrix(c(alpha,as.vector(mubeta),sigma,gam))

  # initilization
  penloglik <-  vector("double", 3)
  loglik <-  vector("double", 3)
  coefficient <- vector("double", length(b))

  penloglik[1] <- out$penloglikset[[index]]
  loglik[1]    <- out$loglikset[[index]]
  for (k in 2:3) {
    ninits <- 1
    maxit <- 2
    # Two EM steps
    out <- cppRegmixPMLE(b, y, x, ztilde, mu0h, sigma0h, m1, p, an, maxit, ninits, epsilon, tau, h, k)
    alpha <- b[1:m1,1] # b0 has been updated
    mubeta <- matrix(b[(1+m1):((q1+1)*m1),1],nrow=q1,ncol=m1)
    sigma <- b[(1+(q1+1)*m1):((q1+2)*m1),1]
    if (!is.null(z)) {
      gam     <- b[((q1+2)*m1+1):((q1+2)*m1+p),1]
    }
    loglik[k]    <- out$loglik[[1]]
    penloglik[k]   <- out$penloglik[[1]]

    # Check singularity: if singular, break from the loop
    if ( any(sigma < 1e-06) || any(alpha < 1e-06) || is.na(sum(alpha)) ) {
      for (i in k:3) {
        loglik[k]    <- -Inf
        penloglik[k]   <- -Inf
      }
      break
    }
    mu.order  <- order(mubeta[1,])
    alpha     <- alpha[mu.order]
    mubeta    <- mubeta[ ,mu.order]
    sigma     <- sigma[mu.order]
    sigma0h <- sigma0h[mu.order]
  }
  coefficient <- as.matrix(c(alpha,as.vector(mubeta),sigma,gam)) # at k=3

  return (list(coefficient = coefficient, loglik = loglik, penloglik = penloglik))
}

#' @description Generates lists of parameters for initial candidates used by
#' the modified EM test for mixture of normals.
#' @export
#' @title regmixPhiInit
#' @name regmixPhiInit
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gam
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gam_1, ..., gam_m))
#' @param z n by p matrix of regressor associated with gam
#' @param h h used as index for pivoting
#' @param tau Tau used to split the h-th component
#' @param ninits number of initial values to be generated
#' @return A list with the following items:
#' \item{alpha}{m+1 by ninits matrix for alpha}
#' \item{mubeta}{q+1 by m+1 by ninits array for mu and beta}
#' \item{sigma}{m+1 by ninits matrix for sigma}
#' \item{gam}{m+1 by ninits matrix for gam}
regmixPhiInit <- function (y, x, z = NULL, parlist, h, tau, ninits = 1)
{
  if (normalregMixtest.env$normalregMix.test.on) # initial values controlled by normalregMix.test.on
    set.seed(normalregMixtest.env$normalregMix.test.seed)

  y  <- as.vector(y)
  n   <- length(y)
  x   <- matrix(x,nrow=n)
  x1   <- cbind(1,x)
  q1   <- ncol(x1)

  alpha0   <- parlist$alpha
  mubeta0  <- parlist$mubeta
  sigma0   <- parlist$sigma
  mu0    <- mubeta0[1,]
  beta0  <- mubeta0[-1,,drop=FALSE]
  gam  <- NULL

  if (!is.null(z)) {
    gam0  <- parlist$gam
    p     <- ncol(z)
    y     <- as.vector(y - z %*% gam0)
    gam1   <- c(rep(1,p),runif(p*(ninits-1),min=-2,max=2))*gam0
    gam   <- matrix(gam1,nrow=p)
  }

  minMU  <- min(y - x%*%beta0)
  maxMU  <- max(y - x%*%beta0)

  m     <- length(alpha0)

  if (m>=2) {
    mid <- (mu0[1:(m-1)]+mu0[2:m])/2  # m-1 by 1
    lb0 <- c(minMU,mid)          # m by 1
    lb  <- c(lb0[1:h],lb0[h:m])      # m+1 by 1
    ub0 <- c(mid,maxMU)          # m by 1
    ub  <- c(ub0[1:h],ub0[h:m])      # m+1 by 1
    beta.hyp <- cbind(beta0[,1:h,drop=FALSE],beta0[,h:m,drop=FALSE])    # q1-1 by m+1
    mu.hyp <- c(mu0[1:h],mu0[h:m])  # m+1 by 1 vector
  } else {  # m=1
    lb  <- c(minMU,minMU)
    ub  <- c(maxMU,maxMU)
    beta.hyp <- cbind(beta0,beta0)   # q1-1 by 2 matirx
    mu.hyp <- c(mu0,mu0)    # 2 by 1 vector
  }

  mubeta <- matrix(0, nrow=q1*(m+1), ncol=ninits)
  for (j in 1:(m+1)) {
    mubeta[(q1*(j-1)+1), ] <- runif(ninits, min=lb, max=ub)
    mubeta[(q1*(j-1)+1), 1] <- mu.hyp[j]
    for (i in 2:q1) {
      mubeta[(q1*(j-1)+i), ] <- beta.hyp[i]*runif(ninits, min=-2, max=2)
      mubeta[(q1*(j-1)+i), 1] <- beta.hyp[i]
    }
  }

  sigma.hyp <- c(sigma0[1:h],sigma0[h:m])  # m+1 by 1
  sigma1 <- c(rep(1,m+1),runif((m+1)*(ninits-1),min=0.25,max=2))*sigma.hyp
  sigma <- matrix(sigma1,nrow=m+1)

  alpha.hyp <- c(alpha0[1:h],alpha0[h:m])  # m+1 by 1
  alpha.hyp[h:(h+1)] <- c(alpha.hyp[h]*tau,alpha.hyp[h+1]*(1-tau))
  alpha <- matrix(rep.int(alpha.hyp,ninits),nrow=m+1)

  list(alpha = alpha, mubeta = mubeta, sigma = sigma, gam = gam)

}  # end function regmixPhiInit

#' @description Estimates parameters of a finite mixture of univariate normals by 
#' the method of penalized maximum likelhood. Using this function is equivalent to
#' calling normalmixPMLE with regressors specified by x as a parameter.
#' @export
#' @title regmixPMLE
#' @name regmixPMLE
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x
#' @param m number of components in the mixture
#' @param z n by p matrix of regressor associated with gam
#' @param vcov.method Method used to compute the variance-covariance matrix, one of \code{"Hessian"} and \code{"OPG"}.
#' The default option is \code{"Hessian"}. When \code{method = "Hessian"}, the variance-covarince matrix is
#' estimated by the Hessian using the formula given in Boldea and Magnus (2009).
#' When \code{method = "OPG"}, the outer product of gradients is used.
#' @param ninits The number of randomly drawn initial values.
#' @param epsilon The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}.
#' @param maxit The maximum number of iterations.
#' @param epsilon.short The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param maxit.short The maximum number of iterations in short EM.
#' @param binit The initial value of parameter vector that is included as a candidate parameter vector
#' @return  A list of class \code{normalMix} with items:
#' \item{coefficients}{A vector of parameter estimates. Ordered as \eqn{\alpha_1,\ldots,\alpha_m,\mu_1,\ldots,\mu_m,\sigma_1,\ldots,\sigma_m,\gam}.}
#' \item{parlist}{The parameter estimates as a list containing alpha, mu, and sigma (and gam if z is included in the model).}
#' \item{vcov}{The estimated variance-covariance matrix.}
#' \item{loglik}{The maximized value of the log-likelihood.}
#' \item{penloglik}{The maximized value of the penalized log-likelihood.}
#' \item{aic}{Akaike Information Criterion of the fitted model.}
#' \item{bic}{Bayesian Information Criterion of the fitted model.}
#' \item{postprobs}{n by m matrix of posterior probabilities for observations}
#' \item{components}{n by 1 vector of integers that indicates the indices of components
#' each observation belongs to based on computed posterior probabilities}
#' \item{call}{The matched call.}
#' \item{m}{The number of components in the mixture.}
#' @note \code{regmixPMLE} maximizes the penalized log-likelihood function
#' using the EM algorithm with combining short and long runs of EM steps as in Biernacki et al. (2003).
#' \code{regmixPMLE} first runs the EM algorithm from \code{ninits}\eqn{* 4m(1 + p)} initial values
#' with the convertence criterion \code{epsilon.short} and \code{maxit.short}.
#' Then, \code{regmixPMLE} uses \code{ninits} best initial values to run the EM algorithm
#' with the convertence criterion \code{epsilon} and \code{maxit}.
#' @references     Biernacki, C., Celeux, G. and Govaert, G. (2003)
#' Choosing Starting Values for the EM Algorithm for Getting the
#' Highest Likelihood in Multivariate Gaussian Mixture Models,
#' \emph{Computational Statistics and Data Analysis}, \bold{41}, 561--575.
#'
#' Boldea, O. and Magnus, J. R. (2009)
#' Maximum Likelihood Estimation of the Multivariate Normal Mixture Model,
#' \emph{Journal of the American Statistical Association},
#' \bold{104}, 1539--1549.
#'
#' Chen, J., Tan, X. and Zhang, R. (2008)
#' Inference for Normal Mixtures in Mean and Variance,
#' \emph{Statistica Sinica}, \bold{18}, 443--465.
#'
#' McLachlan, G. J. and Peel, D. (2000) \emph{Finite Mixture Models}, John Wiley \& Sons, Inc.
#' @examples
#' data(faithful)
#' attach(faithful)
#' regmixPMLE(y = eruptions, x = waiting, m = 1)
#' regmixPMLE(y = eruptions, x = waiting, m = 2)
regmixPMLE <- function (y, x, m = 2, z = NULL, vcov.method = c("Hessian", "OPG", "none"),
                        ninits = 10, epsilon = 1e-08, maxit = 2000,
                        epsilon.short = 1e-02, maxit.short = 500, binit = NULL) {
  y   <- as.vector(y)
  x   <- as.matrix(x)   # n by (q1-1) matrix
  n   <- length(y)
  if (nrow(x) != n) { stop("y and x must have the same number of rows.") }
  x1  <- cbind(1, x)
  q1   <- ncol(x1)

  p       <- 0
  gam   <- NULL
  ninits.short <- ninits*10*(q1+p)*m
  vcov.method <- match.arg(vcov.method)

  if (!is.null(z)) {
    z <- as.matrix(z)
    p <- ncol(z)
    if (nrow(z) != n) { stop("y and z must have the same number of rows.") }
  }

  npar    <- m-1 + (q1+1)*m + p  # number of parameters
  xz      <- cbind(x, z)
  ls.out  <- lsfit(xz, y)
  sd0     <- sqrt(mean(ls.out$residuals^2))

  if (m == 1) {
    mubeta <- as.matrix(unname(ls.out$coeff[1:q1]))
    if (!is.null(z)) {gam <- unname(ls.out$coeff[(q1+1):(q1+p)])}
    res     <- ls.out$residuals
    sigma   <- sqrt(mean(res*res))
    loglik  <- - (n/2)*(1 + log(2*pi) + 2*log(sigma))
    aic     <- -2*loglik + 2*npar
    bic     <- -2*loglik + log(n)*npar
    penloglik <- loglik

    parlist <- list(alpha = 1, mubeta = mubeta, sigma = sigma, gam = gam)
    coefficients <- c(alpha = 1, mubeta = mubeta, sigma = sigma, gam = gam)
    postprobs <- rep(1, n)

  } else {  # m >= 2

    # generate initial values
    tmp <- regmixPMLEinit(y = y, x = x, z = z, ninits = ninits.short, m = m)

    h       <- 0    # setting h=0 gives PMLE
    tau     <- 0.5  # setting tau=0.5 gives PMLE
    k <- 0 # setting k=0 gives PMLE

    sigma0  <- rep(sd0, m)
    mu0     <- double(m)    # dummy
    an      <- 1/n  # penalty term for variance
	
	if (is.null(z))
      ztilde <- matrix(0) # dummy
	else
      ztilde <- z
    
    # short EM
    b0 <- rbind( tmp$alpha, tmp$mubeta, tmp$sigma, tmp$gam )
    if (!is.null(binit)) {
      b0[ , 1] <- binit
    }
    out.short <- cppRegmixPMLE(b0, y, x, ztilde, mu0, sigma0, m, p, an, maxit.short,
                                  ninits.short, epsilon.short)
    # long EM
    components <- order(out.short$penloglikset, decreasing = TRUE)[1:ninits]
    b1 <- b0[ ,components] # b0 has been updated
    out <- cppRegmixPMLE(b1, y, x, ztilde, mu0, sigma0, m, p, an, maxit, ninits, epsilon)

    index     <- which.max(out$penloglikset)
    alpha <- b1[1:m,index] # b0 has been updated
    mubeta <- matrix(b1[(1+m):((q1+1)*m),index],nrow=q1,ncol=m)
    sigma <- b1[(1+(q1+1)*m):((q1+2)*m),index]
    if (!is.null(z)) {
      gam     <- b1[((q1+2)*m+1):((q1+2)*m+p),index]
    }
    penloglik <- out$penloglikset[index]
    loglik    <- out$loglikset[index]
    postprobs <- matrix(out$post[,index], nrow=n)

    aic <- -2*loglik + 2*npar
    bic <- -2*loglik + log(n)*npar

    mu.order  <- order(mubeta[1,])
    alpha     <- alpha[mu.order]
    mubeta    <- mubeta[,mu.order]
    sigma     <- sigma[mu.order]

    postprobs <- postprobs[, mu.order]
    colnames(postprobs) <- c(paste("comp", ".", 1:m, sep = ""))

    mubeta.name <- matrix(0,nrow = q1, ncol = m)
    mubeta.name[1,] <- paste("mu", 1:m, sep = "")

    if (q1 == 2) {
      mubeta.name[2,] <- paste("beta", 1:m,  sep = "")
    } else {
      for (i in 1:(q1-1)) {
        for (j in 1:m) {
          # mubeta.name[i+1,j] <- paste("beta", j, i, sep = "")
          mubeta.name[i+1,j] <- paste("beta", i, ".", j, sep = "")
        }
      }
    }

    parlist <- list(alpha = alpha, mubeta = mubeta, sigma = sigma, gam = gam)
    coefficients <- unlist(parlist)
    names(coefficients)[(m+1):((q1+1)*m)] <- c(mubeta.name)
  }  # end m >= 2

  if (vcov.method == "none") {
    vcov <- NULL
  } else {
    vcov <- regmixVcov(y = y, x = x, coefficients = coefficients, z = z , vcov.method = vcov.method)
  }

  a <- list(coefficients = coefficients, parlist = parlist, vcov = vcov, loglik = loglik,
            penloglik = penloglik, aic = aic, bic = bic, postprobs = postprobs,
            components = getComponentcomponents(postprobs),
            call = match.call(), m = m, label = "PMLE")

  class(a) <- "normalregMix"

  a

}  # end function regmixPMLE

#' Generate initial values used by the PMLE of mixture of normals
#' @export
#' @title regmixPMLEinit
#' @name regmixPMLEinit
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x
#' @param z n by p matrix of regressor associated with gam
#' @param ninits number of initial values to be generated
#' @param m The number of components in the mixture
#' @return A list with the following items:
#' \item{alpha}{m by ninits matrix for alpha}
#' \item{mubeta}{q+1 by m by ninits array for mu and beta}
#' \item{sigma}{m by ninits matrix for sigma}
#' \item{gam}{m by ninits matrix for gam}
regmixPMLEinit <- function (y, x, z = NULL, ninits = 1, m = 2)
{
  if (normalregMixtest.env$normalregMix.test.on) # initial values controlled by normalregMix.test.on
    set.seed(normalregMixtest.env$normalregMix.test.seed)

  n  <- length(y)
  q1  <- ncol(x)+1
  p  <- ncol(z)

  gam <- NULL
  if (!is.null(z)) {
    out     <- lsfit(cbind(x, z), y)
    gam0  <- out$coef[(q1+1):(q1+p)]
    gam   <- matrix(runif(p*ninits, min=0.5, max=1.5), nrow=p)*gam0
    mubeta_hat <- out$coef[1:q1]
    y     <- y - z %*% gam0
    r     <- out$residuals
    stdR  <- sd(r)
  } else {
    out         <- lsfit(x, y)
    mubeta_hat  <- out$coef
    r           <- out$residuals
    stdR        <- sd(r)
  }

  alpha <- matrix(runif(m*ninits), nrow=m)
  alpha <- t(t(alpha)/colSums(alpha))


  minMU <- min(y - x %*% mubeta_hat[-1])
  maxMU <- max(y - x %*% mubeta_hat[-1])
  mubeta <- matrix(0, nrow=q1*m, ncol=ninits)
  for (j in 1:m) {
    mubeta[(q1*(j-1)+1), ] <- runif(ninits, min=minMU, max=maxMU)
    for (i in 2:q1) {
      mubeta[(q1*(j-1)+i), ] <- mubeta_hat[i]*runif(ninits, min=-2, max=2)
    }
  }
  sigma <- matrix(runif(m*ninits, min=0.01, max=1), nrow=m)*stdR

  list(alpha = alpha, mubeta = mubeta, sigma = sigma, gam = gam)

}  # end function regmixPMLEinit

#' Computes the variance-covariance matrix of the MLE of m-component normal mixture.
#' @export
#' @title regmixVcov
#' @name regmixVcov
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x
#' @param coefficients (alpha_1, ..., alpha_m, mu_1, ..., mu_m, sigma_1, ..., sigma_m, gam)
#' @param z n by p matrix of regressor associated with gam
#' @param vcov.method Method used to compute the variance-covariance matrix,
#' one of \code{"Hessian"} and \code{"OPG"}. #' The default option is \code{"Hessian"}.
#' When \code{method = "Hessian"}, the variance-covarince matrix is
#' estimated by the Hessian using the formula given in Boldea and Magnus (2009).
#' When \code{method = "OPG"}, the outer product of gradients is used.
#' @return The variance-covariance matrix of the MLE of
#' m-component normal mixture given the data and coefficients.
#' @references   Boldea, O. and Magnus, J. R. (2009)
#' Maximum Likelihood Estimation of the Multivariate Normal Mixture Model,
#' \emph{Journal of the American Statistical Association},
#' \bold{104}, 1539--1549.
regmixVcov <- function(y, x, coefficients, z = NULL, vcov.method = c("Hessian", "OPG")) {
  # Computes the variance-covariance matrix of the MLE of m-component normal regression mixture
  # Input
  #  y  : n by 1 vector of dependent variable
  #  x  : n by q1-1 matrix of regressor NOT including an intercept
  #  coefficients : (alpha_1,...,alpha_m,mubeta_1^T, ...,mubeta_m^T,sigma_1, ..., sigma_m,gam^T)
  #  z  : n by p matrix of regressor associated with gam
  # Output
  #  vcov: variance-covariance matrix

  y     <- as.vector(y)
  n     <- length(y)
  len   <- length(coefficients)
  p     <- 0
  gam <- NULL
  vcov.method <- match.arg(vcov.method)

  if (!is.null(z)) {
    z <- as.matrix(z)
    p <- ncol(z)
    gam <- coefficients[(len-p+1):len]
  }

  x   <- as.matrix(x)
  x1  <- cbind(1, x)
  q1   <- ncol(x1)
  q   <- ncol(x)

  m  <- (len-p)/(3+q)
  if (round(m) != m)
    stop("The dimension of the coefficients is incompatible with x and z. Please check the data.")

  alpha   <- coefficients[1:m]  # m-vector
  mubeta  <- matrix(coefficients[(m+1):((2+q)*m)], nrow=q+1, ncol=m)  # q+1 by m
  sigma   <- coefficients[((2+q)*m+1):((3+q)*m)]  # m-vector

  if (m == 1) {
    xz1 <- cbind(x1,z)
    I <- matrix(0, nrow=q1+p+1, ncol=q1+p+1)
    I[1:(q1+p), 1:(q1+p)] <- t(xz1) %*% xz1/sigma^2
    I[(q1+p+1), (q1+p+1)] <- n*2/sigma^2
    if (p != 0){
      s.1 <- c(1:q1, (q1+p+1), (q1+1):(q1+p))
      I <- I[s.1, s.1]
    }

    vcov <- solve(I)
    # Because the variance is parameterized as sigma^2, we convert it to sigma
    c.mat.vec <- c(rep(1,q1),(1/sigma^(1/2))/2,rep(1,p))
    vcov <- diag(c.mat.vec) %*% vcov %*% diag(c.mat.vec)

  } else {  # end of if (m == 1)
    # m >= 2
    # Compute posterior probabilities, and adjust y if z is present
    sigma0  <- rep(1, m)  # dummy
    mu0     <- double(m)  # dummy
    an      <- 1/n  # penalty term for variance
    h       <- 0
    tau     <- 0.5
    k <- 0
    epsilon <- 1e-08
	maxit = 2
    ninits = 1
    b <- matrix( rep( coefficients, ninits), ncol = ninits)
	
    if (is.null(z)) 
      out.p <- cppRegmixPMLE(b, y, x, matrix(0), mu0, sigma0, m, p, an, maxit, ninits, epsilon, tau, h, k)
    else 
	  {
      out.p <- cppRegmixPMLE(b, y, x, z, mu0, sigma0, m, p, an, maxit, ninits, epsilon, tau, h, k)
      # Adjust y
      y <- as.vector(y - z %*% gam)
    }

    post <- matrix(out.p$post, nrow=n)

    p2 <- seq(q1+1, (q1+1)*m, by=q1+1)  # sequence of q1+1, (q1+1)*2, ... , (q1+1)*m
    p1 <- (1:((q1+1)*m))[-p2]        # other values from 1, ..., (q1+1)*m

    a <- diag(1/alpha[-m], nrow=m-1, ncol=m-1)
    a <- cbind(a, -1/alpha[m])  # m-1 by m matrix of a_i's
    abar <- a %*% t(post)  # m-1 by n

    xtheta <- x1 %*% mubeta  # n by m

    Z0 <- t(t(y-xtheta)/sigma)          # normalized data, n by m
    f <- t(t(exp(-Z0^2/2)/sqrt(2*pi))/sigma)  # pdf, n by m
    phi <- t(t(f)*alpha)            # n by m
    f0 <- rowSums(phi)              # data pdf, n by 1

    vinv <- 1/(sigma*sigma)  # m-vector

    b <- t(t(Z0)/sigma)  # n by m
    B <- t(vinv - t(b*b))  # n by m

    c0 <- array(0,dim=c(n, m, q1+1))
    c0[, , (1:q1)] <- array(tKR(x1, b), dim=c(n, m, q1))
    c0[, , (q1+1)] <- -B/2

    # Compute Hessian-based I
    if (vcov.method == "Hessian") {
      other.method = "OPG"
      C0 <- array(0, dim=c(n, m, q1+1, q1+1))
      x11 <- array(tKR(x1, x1), dim = c(n, q1, q1))
      for (i in 1:m) {
        C0[, i, (1:q1), (1:q1)] <- x11*vinv[i]
        C0[, i, (1:q1), q1+1]   <- C0[, i, q1+1, (1:q1)] <- x1*b[, i]*vinv[i] # n by q1
      }
      C0[, , q1+1, q1+1] <- t((vinv - 2*t(B))*vinv)/2      # n by m

      Q.pi <- - abar %*% t(abar)  # m-1 by m-1

      Q.pi.theta <- matrix(0,nrow=m-1,ncol=(q1+1)*m)  # m-1 by (q1+1)*m
      for (i in 1:m) {
        zi <- a[, i] - abar  # m-1 by n
        wi <- c0[, i, ]*post[, i]  # n by q1+1
        Q.i <- colSums(tKR(wi, t(zi)))  # (q1+1)*(m-1) vector
        # first q1*(m-1) elements correspond to mubeta x pi,
        # last m-1 elements correspond to sigma x pi,
        Q.pi.theta[,(q1*(i-1)+1):(q1*i)] <- matrix(Q.i[1:(q1*(m-1))],ncol=q1)  # m-1 by q1 matrix
        Q.pi.theta[, q1*m+i] <- Q.i[(q1*(m-1)+1):((q1+1)*(m-1))]  # m-1 vector
      }

      Q.theta <- matrix(0, nrow=(q1+1)*m, ncol=(q1+1)*m)
      for (i in 2:m) {  # off-diagonal blocks
        for (j in 1:(i-1)) {
          wi  <- c0[, i, ]*post[, i] # n by q1+1
          wj  <- c0[, j, ]*post[, j] # n by q1+1
          Q.ij <- - colSums(tKR(wi, wj))  # (q1+1)*(q1+1) vector
          Q.theta[((q1+1)*(i-1)+1):((q1+1)*i), ((q1+1)*(j-1)+1):((q1+1)*j)] = t(matrix(Q.ij, nrow=q1+1, ncol=q1+1))
        }
      }

      Q.theta <- Q.theta + t(Q.theta)
      for (i in 1:m) {  # diagonal blocks
        C.ii   <- array(C0[, i, , ], dim=c(n, q1+1, q1+1))
        Q.ii.1   <- apply(C.ii*post[,i], c(2, 3), sum)
        w.ii   <- tKR(c0[, i, ], c0[, i, ])*post[, i]*(1-post[, i])
        Q.ii.2   <- matrix(colSums(w.ii), nrow=q1+1, ncol=q1+1)
        Q.theta[((q1+1)*(i-1)+1):((q1+1)*i), ((q1+1)*(i-1)+1):((q1+1)*i)] <- -Q.ii.1 + Q.ii.2
      }

      # q1+1,2*(q1+1),...,m*(q1+1)th rows and columns = sigma
      # other rows and columns = mubeta
      Q.theta <- Q.theta[c(p1, p2), c(p1, p2)]  # first block = wrt mubeta, second blosk = wrt sigma

      dimI <- m-1+(q1+1)*m
      I <- matrix(0, nrow=dimI, ncol=dimI)
      I[1:(m-1), 1:(m-1)] <- - Q.pi
      I[1:(m-1), m:dimI]  <- - Q.pi.theta
      I[m:dimI, 1:(m-1)]  <- - t(Q.pi.theta)
      I[m:dimI, m:dimI]   <- - Q.theta

      if (!is.null(z)) {
        dbar <-  z*rowSums(post*b)  # n by p
        Q.gam.theta <- matrix(0, nrow=p, ncol=(q1+1)*m)  # p by (q1+1)*m matrix
        for (i in 1:m) {
          C.i <- array(C0[, i, 1, ], dim=c(n, q1+1))  # n by q1+1
          Q.i.1 <- colSums(tKR(-C.i+b[, i]*c0[, i, ], z*post[, i])) # p*(q1+1) vector
          Q.i.2 <- colSums(tKR(c0[, i, ]*post[, i], dbar))  # p*(q1+1) vector
          Q.gam.theta[, ((q1+1)*(i-1)+1):((q1+1)*i)] <- matrix(Q.i.1+Q.i.2, nrow=p, ncol=q1+1)
        }

        Q.gam.theta <- Q.gam.theta[, c(p1, p2), drop=FALSE]  # p by (q1+1)*m
        w1 <- (post*b)%*%t(a) - rowSums(post*b)*t(abar)  # n by m-1
        Q.pi.gam.0 <- colSums(tKR(w1, z))  # (m-1)*p vector
        Q.pi.gam  <- matrix(Q.pi.gam.0, nrow=m-1, ncol=p)
        Q.gam     <- - t(z)%*%(z*rowSums(post*B)) -
          matrix(colSums(tKR(dbar, dbar)), nrow=p, ncol=p)

        I <- cbind(I, -rbind(Q.pi.gam, t(Q.gam.theta)))
        I <- rbind(I, -cbind(t(Q.pi.gam), Q.gam.theta, Q.gam))
      }  # end if (!is.null(z))

    }  else {  # compute I with (method == "OPG")
      other.method = "Hessian"
      c0.a <- array(0, dim=c(n, m, 2))
      c0.a[, , 1] <- b  # n by m
      c0.a[, , 2] <- -B/2  # n by m

      score <- t(abar)

      for (j in 1:m) {
        # score.o <- cbind(score.o, c0[, j, ]*post[, j])
        score <- cbind(score, x1*c0.a[, j, 1]*post[, j], c0.a[, j, 2]*post[, j])
        # print(all.equal(score.o, score))
      }

      ind <- c(1:(m-1), p1+m-1, p2+m-1)
      score <- score[, ind]
      I <- t(score) %*% score

      if (!is.null(z))  {
        dbar <-  z*rowSums(post*b)  # n by p
        score <- cbind(score, dbar)
        I <- t(score) %*% score
      }

    }  # end if (method=="OPG")

    vcov <- try(solve(I))
    if (class(vcov) == "try-error" || any(diag(vcov) <0) ) {
      vcov <- matrix(NaN, nrow = (2+q1)*m-1+p, ncol = (2+q1)*m-1+p)
      warning("Fisher information matrix is singular and/or the
              variance is estimated to be negative. Consider using vcov.method=\"",other.method,"\".")
    }

    # Because the variance is parameterized as sigma^2, we convert it to sigma

    c.mat.vec <- c(rep(1, m-1+m*q1), (1/sigma^(1/2))/2, rep(1, p))
    vcov <- diag(c.mat.vec) %*% vcov %*% diag(c.mat.vec)
    # vcov.opg <- diag(c.mat.vec) %*% vcov.opg %*% diag(c.mat.vec)

    # Add the variance of alpha_m
    M.mat <- diag(len-1)
    M.mat <- rbind(M.mat[1:(m-1),], c(rep(-1,m-1),rep(0,len-m)), M.mat[m:(len-1),])

    vcov <- M.mat %*% vcov %*% t(M.mat)
    # vcov.opg <- M.mat %*% vcov.opg %*% t(M.mat)

  }   # end else (i.e., m >= 2)

  vcov

}  # end function regmixVcov


#' @description Computes residuals used in function LR_2.comp.
#' @export
#' @title obj_zIz
#' @name obj_zIz
#' @param b n by m matrix of normalized data.
#' @param Z nrep by dim(\eqn{q_\lambda}) matrix of random vectors.
#' @param I dim(\eqn{q_\lambda}) by dim(\eqn{q_\lambda}) matrix of the the inverse of the covariance matrix of Z.
#' @return nrep by dim(\eqn{q_\lambda}) matrix of residuals.
#' @keywords{intenal}
obj_zIz <- function(b,Z,I) {
  q       <- length(b)-2
  lam_mu  <- b[1]
  lam_sig <- b[2]
  lam_bet <- b[3:length(b)]
  t       <- double(2+2*q+q*(q+1)/2)
  t[1]          <- 6*lam_mu*lam_sig
  t[2:(q+1)]    <- 2*lam_bet*lam_mu
  t[q+2]        <- 12*lam_sig^2
  t[(q+3):(2+2*q)]  <- 2*lam_bet*lam_sig

  t[(2+2*q+1):(2+2*q+q)] <- lam_bet^2

  if (q>=2) {
    ii <- combn(1:q,2)
    t[(3+3*q):(2+2*q+q*(q+1)/2)] <- lam_bet[ii[1,]]*lam_bet[ii[2,]]
  }

  R <- chol(I)
  res <- R %*% (t-as.vector(Z))

  res

} # end function obj_zIz

#' @description Computes Jacobian of the residuals used function LR_2.comp.
#' @export
#' @title obj_zIz.jac
#' @name obj_zIz.jac
#' @param b dim(\eqn{\lambda}) by 1 vector of \eqn{\lambda}.
#' @param Z nrep by dim(\eqn{q_\lambda}) matrix of random vectors.
#' @param I dim(\eqn{q_\lambda}) by dim(\eqn{q_\lambda}) matrix of the the inverse of the covariance matrix of Z.
#' @return dim(\eqn{q_\lambda}) by dim(\eqn{q_\lambda}) matrix of Jacobian.
#' @keywords{intenal}
obj_zIz.jac <- function(b,Z,I) {

  q     <- length(b)-2
  lam_mu  <- b[1]
  lam_sig <- b[2]
  lam_bet <- b[3:length(b)]
  t     <- double(2+2*q+q*(q+1)/2)
  t[1]        <- 6*lam_mu*lam_sig
  t[2:(q+1)]  <- 2*lam_bet*lam_mu
  t[q+2]      <- 12*lam_sig^2
  t[(q+3):(2+2*q)]  <- 2*lam_bet*lam_sig

  t[(2+2*q+1):(2+2*q+q)]  <- lam_bet^2

  if (q >= 2) {
    ii <- combn(1:q,2)
    t[(3+3*q):(2+2*q+q*(q+1)/2)] <- lam_bet[ii[1,]]*lam_bet[ii[2,]]
  }

  R <- chol(I)

  R_sub = NULL
  if (q >= 2) {
    R_sub = cbind(lam_bet[2:q], lam_bet[1]*diag(q-1))
    if (q >= 3) {
      for (i in 2:(q-1)) {
        R_sub_2 = cbind(matrix(0,nrow=q-i,ncol=i-1),lam_bet[(i+1):q],lam_bet[i]*diag(q-i))
        R_sub = rbind(R_sub, R_sub_2)
      }

    }

  }

  G <- cbind(c(6*lam_sig, 2*lam_bet, double(1+q+q*(q+1)/2)),
             c(6*lam_mu, double(q), 24*lam_sig, 2*lam_bet, double(q*(q+1)/2)),
             rbind(double(q), 2*lam_mu*diag(q), double(q), 2*lam_sig*diag(q),
                   2*diag(lam_bet,nrow=length(lam_bet)), R_sub))

  jac <- R%*%G

  return(jac)

} # end function obj_zIz.jac

#' @description Computes LR_2 in computing the critcal values of regmixMEMtest.
#' @export
#' @title LR_2.comp
#' @name LR_2.comp
#' @param Z nrep by dim(\eqn{\lambda}) matrix of random vectors.
#' @param I dim(\eqn{\lambda}) by dim(\eqn{\lambda}) matrix of the the inverse of the covariance matrix of Z.
#' @param q Number of elements in x.
#' @param ninits Number of initial values used in maximizing the objective function
#' @return nrep by 1 vector of the values of the test statisic.
#' @keywords{intenal}
LR_2.comp <- function(Z, I, q, ninits = 25) {
  if (normalregMixtest.env$normalregMix.test.on) # initial values controlled by normalregMix.test.on
    set.seed(normalregMixtest.env$normalregMix.test.seed)

  con.nls = minpack.lm::nls.lm.control(maxit = 400)

  n_t <- q+2
  ll <- vector("list",n_t)
  for (i in 1:(q+2)) {ll[[i]] <- c(1,-1)}
  # g <- t(as.matrix(expand.grid(ll)))
  g <- t((expand.grid(ll)))

  lam_sig_0 <- sqrt(abs(Z[2+q]/12))
  lam_bet_0 <- sqrt(abs(Z[(2+2*q+1):(2+2*q+q)]))
  lam_mu_0 <- sqrt(abs(Z[1]*Z[2]/12/lam_sig_0/lam_bet_0[1]))
  t0 <- cbind(c(lam_mu_0,lam_sig_0,lam_bet_0), 8*matrix(rnorm(n_t*(ninits-1)), nrow=n_t))

  LR_2_all <- double(ninits)
  LR_2_all.optim <- double(ninits)

  for (irep in 1:ninits) {
    tg <- g*t0[,irep]  # q+2 by ncol(g) matrix
    res <- apply(tg,2,obj_zIz,Z,I)
    obj_val <- colSums(res^2)
    lm.par <- tg[,which.min(colSums(res^2))]  # uses the column of tg that gives the smallest obj value
    nls.out <- minpack.lm::nls.lm(par=lm.par, fn = obj_zIz, jac = obj_zIz.jac, control = con.nls, Z=Z, I=I)
    LR_2_all[irep] <- sum((nls.out$fvec)^2)
  }

  LR_2 <- Z %*% I %*% Z - min(LR_2_all)

  LR_2

} # end function LR_2.comp
