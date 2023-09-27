#' @description Computes the variance-covariance matrix of the MLE of 
#' m-component normal mixture.
#' @export
#' @title normalmixVcov
#' @name normalmixVcov
#' @param y n by 1 vector of data.
#' @param coefficients parameter vector ordered as (alpha_1, ..., alpha_m, mu_1, ..., mu_m, sigma_1, ..., sigma_m, gamma).
#' @param z n by p matrix of regressor associated with gamma.
#' @param vcov.method method used to compute the variance-covariance matrix,
#' one of \code{"Hessian"} and \code{"OPG"}. Default is \code{"Hessian"}.
#' When \code{method = "Hessian"}, the variance-covarince matrix is
#' estimated by the Hessian using the formula given in Boldea and Magnus (2009).
#' When \code{method = "OPG"}, the outer product of gradients is used.
#' @return variance-covariance matrix of the MLE of
#' m-component normal mixture given the data and coefficients.
#' @references  Boldea, O. and Magnus, J. R. (2009)
#' Maximum Likelihood Estimation of the Multivariate Normal Mixture Model,
#' \emph{Journal of the American Statistical Association},
#' \bold{104}, 1539--1549.
normalmixVcov <- function(y, coefficients, z = NULL, vcov.method = c("Hessian", "OPG"))
{
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
      I <- n*diag(c(1/sigma^2, 0.5/sigma^4))  # information matrix when
      # the variance is parameterized as sigma^2
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
    c.mat.vec <- c(1, 0.5/sigma, rep(1, p))
    vcov <- diag(c.mat.vec) %*% vcov %*% diag(c.mat.vec)

  } else { # end of if (m == 1)
    # m >= 2

    # Compute posterior probabilities, and adjust y if z is present
    sigma0  <- rep(1, m)  # dummy
    mu0     <- double(m)  # dummy
    an      <- 1/n  # penalty term for variance
    cn      <- 0   # penalty term for tau
    h       <- 0
    tau     <- 0.5
    k       <- 0
    epsilon <- 1e-08
    maxit = 1
    ninits = 1

  # b <- matrix( rep( coefficients, ninits), ncol = ninits)
  # if (is.null(z)) {
  #   out.p <- cppNormalmixPMLE(b, y, mu0, sigma0, m, an, cn, maxit, ninits, epsilon, tau, h, k)
  # } else {
  # out.p <- cppNormalmixPMLE(b, y, z,mu0, sigma0, m, p, an, cn, maxit, ninits, epsilon, tau, h, k)
  # # Adjust y
  # y <- y - z %*% gamma
  # }
  # post <- matrix(out.p$post, nrow=n)
  
    if (!is.null(z)) {
      y <- y - z %*% gamma  # Adjust y
    }

    post <- matrix(0, nrow=n, ncol=m)
    for (j in 1:m){
      post[,j] <- alpha[j]*dnorm(y, mean = mu[j], sd = sigma[j])
    }
    f <- post %*% rep(1,m)
    post <- post/(f %*% rep(1,m))
    
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
    if (any(class(vcov) == "try-error") || any(diag(vcov) <0) ) {
      vcov <- matrix(NaN, nrow = 3*m-1+p, ncol = 3*m-1+p)
      warning("Fisher information matrix is singular and/or the
              variance is estimated to be negative. Consider using vcov.method=\"",other.method,"\".")
    }

    # Because the variance is parameterized as sigma^2, we convert is to sigma

    c.mat.vec <- c(rep(1, m-1+m), 0.5/sigma, rep(1, p))
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


#' @description Estimates parameters of a finite mixture of univariate normals by 
#' the penalized maximum likelhood estimator.
#' @export
#' @title normalmixPMLE
#' @name normalmixPMLE
#' @param y n by 1 vector of data.
#' @param m number of components in the mixture.
#' @param z n by p matrix of regressor associated with gamma.
#' @param vcov.method method used to compute the variance-covariance matrix, one of \code{"Hessian"} and \code{"OPG"}.
#' Default is \code{"Hessian"}. When \code{method = "Hessian"}, the variance-covarince matrix is
#' estimated by the Hessian using the formula given in Boldea and Magnus (2009).
#' When \code{method = "OPG"}, the outer product of gradients is used.
#' @param ninits number of randomly drawn initial values.
#' @param epsilon convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}.
#' @param maxit maximum number of iterations.
#' @param epsilon.short convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param maxit.short maximum number of iterations in short EM.
#' @param binit initial value of parameter vector that is included as a candidate parameter vector.
#' @return  A list of class \code{normalregMix} with items:
#' \item{coefficients}{vector of parameter estimates. Ordered as \eqn{\alpha_1,\ldots,\alpha_m,\mu_1,\ldots,\mu_m,\sigma_1,\ldots,\sigma_m,\gamma}.}
#' \item{parlist}{parameter estimates as a list containing alpha, mu, and sigma (and gam if z is included in the model).}
#' \item{vcov}{estimated variance-covariance matrix.}
#' \item{loglik}{maximized value of the log-likelihood.}
#' \item{penloglik}{maximized value of the penalized log-likelihood.}
#' \item{aic}{Akaike Information Criterion of the fitted model.}
#' \item{bic}{Bayesian Information Criterion of the fitted model.}
#' \item{postprobs}{n by m matrix of posterior probabilities for observations.}
#' \item{components}{n by 1 vector of integers that indicates the indices of components
#' each observation belongs to based on computed posterior probabilities.}
#' \item{call}{The matched call.}
#' \item{m}{number of components in the mixture.}
#' @note \code{normalmixPMLE} maximizes the penalized log-likelihood function
#' (Chen et al., 2008) 
#' using the EM algorithm with combining short and long runs of EM steps as in Biernacki et al. (2003).
#' \code{normalmixPMLE} first runs the EM algorithm from \code{ninits}\eqn{* 4m(1 + p)} initial values
#' with the convergence criterion \code{epsilon.short} and \code{maxit.short}.
#' Then, \code{normalmixPMLE} uses \code{ninits} best initial values to run the EM algorithm
#' with the convergence criterion \code{epsilon} and \code{maxit}.
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
# McLachlan, G. J. and Peel, D. (2000) \emph{Finite Mixture Models}, John Wiley \& Sons, Inc.
#' @examples
#' data(faithful)
#' attach(faithful)
#'
#' normalmixPMLE(y = waiting, m = 1)
#' normalmixPMLE(y = waiting, m = 2)
#'
#' out <- normalmixPMLE(y = waiting, m = 2)
#' summary(out)
normalmixPMLE <- function (y, m = 2, z = NULL, vcov.method = c("Hessian", "OPG", "none"),
                           ninits = 25, epsilon = 1e-08, maxit = 2000,
                           epsilon.short = 1e-02, maxit.short = 500, binit = NULL) {
  y <- as.vector(y)
  n <- length(y)
  p <- 0
  gamma <- NULL
  ninits.short <- ninits*10*(1 + p)*m
  vcov.method <- match.arg(vcov.method)

  if (!is.null(z)) {
    z     <- as.matrix(z)
    p     <- ncol(z)
    if (nrow(z) != n) { stop("y and z must have the same number of rows.") }
    ls.out   <- lsfit(z, y)
    sd0    <- sqrt(mean(ls.out$residuals^2))
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

    if (is.null(z)) {
      # generate initial values for short EM
      tmp <- normalmixPMLEinit(y = y, z = z, ninits = ninits.short, m = m)

      # the following values for (h, k, tau, an) are given by default
      # h       <- 0  # setting h=0 gives PMLE
      # k       <- 0  # k is set to 0 because this is PMLE
      # tau     <- 0.5  # tau is set to 0.5 because this is PMLE
      an      <- 1/n  # penalty term for variance
      cn      <- 0   # penalty term for tau
      sigma0  <- rep(sd0, m)
      mu0     <- double(m+1) # dummy

      ztilde <- matrix(0) # dummy
      # short EM
      b0 <- as.matrix(rbind( tmp$alpha, tmp$mu, tmp$sigma))
      if (!is.null(binit)) {
        b0[ , 1] <- binit
      }
      out.short <- cppNormalmixPMLE(b0, y, mu0, sigma0, m, an, cn, maxit.short,
                                ninits.short, epsilon.short)
      # long EM
      components <- order(out.short$penloglikset, decreasing = TRUE)[1:ninits]
      b1 <- b0[ ,components] # b0 has been updated
      out <- cppNormalmixPMLE(b1, y, mu0, sigma0, m, an, cn, maxit, ninits, epsilon)
    
      index     <- which.max(out$penloglikset)
      alpha <- b1[1:m,index] # b0 has been updated
      mu <- b1[(1+m):(2*m),index]
      sigma <- b1[(1+2*m):(3*m),index]
      gamma <- NULL
      penloglik <- out$penloglikset[index]
      loglik    <- out$loglikset[index]
      postprobs <- matrix(out$post[,index], nrow=n)
    } else {
      # generate initial values for nlopt/SLSQP
      ninits.nlopt <- ninits*(1 + p*p)*m
      tmp <- normalmixPMLEinit(y = y, z = z, ninits = ninits.nlopt, m = m)
      h       <- 0  # dummy
      tau     <- 0.5
      an      <- 1/n  # penalty term for variance
      cn      <- 0   # penalty term for tau
      sigma0  <- rep(sd0, m)
      mu0     <- double(m+1) # dummy
      
      b0 <- as.matrix(rbind( tmp$alpha, tmp$mu, tmp$sigma, tmp$gamma))
      if (!is.null(binit)) {
        b0[ , 1] <- binit
      }
      out.list <- vector(mode='list', length=ninits)
      opts <- list("algorithm"="NLOPT_LD_SLSQP", "xtol_rel"=1.0e-8, "maxeval"=1000)
      lb <- c(rep(0,m), rep(-Inf,m), rep(0,m), rep(-Inf,p))
      ub <- c(rep(1,m), rep(Inf,m), rep(Inf,m), rep(Inf,p))
      
      suppressWarnings({
        for (jj in 1:ninits.nlopt){
          if (is.finite(pllnormalmix_z_cpp(b0[,jj], y, z, m, p, an, cn, sigma0, tau, h)) &&
              all(is.finite(pllnormalmix_z_grad_cpp(b0[,jj], y, z, m, p, an, cn, sigma0, tau, h))) ) {
              # out.list[[jj]] <- nloptr(x0=b0[,jj], eval_f=pllnormalmix_z,
              out.list[[jj]] <- nloptr(x0=b0[,jj], eval_f=pllnormalmix_z_cpp,
                                        eval_grad_f=pllnormalmix_z_grad_cpp,
                                        lb=lb, ub=ub, 
                                        eval_g_eq = normalmix_heq_cpp,
                                        eval_jac_g_eq = normalmix_heq_grad_cpp,
                                        opts=opts,
                                        y=y, z=z, m=m, p=p, an=an, cn=cn,
                                        sigma0=sigma0, tau=tau, h=h)
          } else {
            out.list[[jj]] <- list(objective = Inf)
          } # enf of if (is.finite...)
        } # end of for (jj in 1:ninits) loop
      })
      
      penloglikset <- sapply(out.list,"[[", "objective")
      index <- which.min(penloglikset)
      out <- out.list[[index]]
      penloglik <- -out$objective
      
      theta <- out$solution
      alpha  <- theta[1:m]
      mu     <- theta[(m+1):(2*m)]
      sigma  <- theta[(2*m+1):(3*m)]
      gamma  <- matrix(theta[(3*m+1):(3*m+p)], ncol=1)
      ss <- (sigma0/sigma)^2
      loglik <- penloglik + an*sum(ss + log(1/ss) -1)

      ytilde <- y - z %*% gamma
      
      postprobs <- matrix(0, nrow=n, ncol=m)
      for (j in 1:m){
        postprobs[,j] <- alpha[j]*dnorm(ytilde, mean = mu[j], sd = sigma[j])
      }
      f <- postprobs %*% rep(1,m)
      postprobs <- postprobs/(f %*% rep(1,m))
    
    }
    
    aic     <- -2*loglik + 2*(m-1+2*m+p)
    bic     <- -2*loglik + log(n)*(m-1+2*m+p)

    mu.order  <- order(mu)
    alpha     <- alpha[mu.order]
    mu        <- mu[mu.order]
    sigma     <- sigma[mu.order]

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
            penloglik = penloglik, aic = aic, bic = bic, 
            postprobs = postprobs, components = getComponentcomponents(postprobs),
            call = match.call(), m = m, label = "PMLE")

  class(a) <- "normalregMix"

  a
}  # end function normalmixPMLE

#' @description Generates lists of parameters for initial candidates used by
#' \code{normalmixMaxPhi}.
#' @title normalmixPhiInit
#' @name normalmixPhiInit
#' @param y n by 1 vector of data.
#' @param parlist parameter estimates as a list containing alpha, mu, sigma,
#' and gamma in the form of (alpha = (alpha_1, ..., alpha_m),
#' mu = (mu_1, ..., mu_m), sigma = (sigma_1, ..., sigma_m), gam).
#' @param z n by p matrix of regressor associated with gamma.
#' @param h h used as index for pivoting.
#' @param tau tau used to split the h-th component.
#' @param ninits number of initial values to be generated.
#' @return A list with the following items:
#' \item{alpha}{m+1 by ninits matrix for alpha.}
#' \item{mu}{m+1 by ninits matrix for mu.}
#' \item{sigma}{m+1 by ninits matrix for sigma.}
#' \item{gamma}{p by ninits matrix for gamma.}
normalmixPhiInit <- function (y, parlist, z = NULL, h, tau, ninits = 1)
{
  if (normalregMixtest.env$normalregMix.test.on) # initial values controlled by normalregMix.test.on
  set.seed(normalregMixtest.env$normalregMix.test.seed)
  
  n     <- length(y)
  p     <- ncol(z)

  mu0      <- parlist$mu
  sigma0   <- parlist$sigma
  alpha0   <- parlist$alpha
  m       <- length(alpha0)

  if (is.null(z))
    gamma <- NULL
  else {
    z <- matrix(z, ncol=p)
    gamma0  <- matrix(parlist$gamma, ncol=1)
    y    <- y - z %*% gamma0
    gamma <- matrix(runif(p*ninits,min=0.5,max=1.5),nrow=p)*(gamma0 %*% rep(1,ninits))
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


#' @description Generate initial values used by \code{normalmixPMLE}.
#' @title normalmixPMLEinit
#' @name normalmixPMLEinit
#' @param y n by 1 vector of data.
#' @param z n by p matrix of regressor associated with gamma.
#' @param ninits number of initial values to be generated.
#' @param m number of components in the mixture.
#' @return A list with the following items:
#' \item{alpha}{m by ninits matrix for alpha.}
#' \item{mu}{m by ninits matrix for mu.}
#' \item{sigma}{m by ninits matrix for sigma.}
#' \item{gamma}{p by ninits matrix for gamma.}
normalmixPMLEinit <- function (y, z = NULL, ninits = 1, m = 2)
{
  if (normalregMixtest.env$normalregMix.test.on) # initial values controlled by normalregMix.test.on
    set.seed(normalregMixtest.env$normalregMix.test.seed)
  
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

  a <- kmeans(y, m)
  alpha1s <- a$size/n
  mu1s <- a$centers
  sigma1s <- sqrt(pmax(a$withinss/a$size, 0.01))
  
  mu1sorder <- order(mu1s)
  alpha1s <- matrix(alpha1s[mu1sorder], nrow=m)
  mu1s <- matrix(mu1s[mu1sorder], nrow=m)
  sigma1s <- matrix(sigma1s[mu1sorder], nrow=m)
  alpha[,1] <- alpha1s
  mu[,1] <- mu1s
  sigma[,1] <- sigma1s
  
  list(alpha = alpha, mu = mu, sigma = sigma, gamma = gamma)

}  # end function normalmixPMLEinit

#' @description Generates mixed normal random variables with regressor x.
#' @export
#' @title rnormregmix
#' @name rnormregmix
#' @param n number of observations.
#' @param alpha m by 1 vector that represents proportions of components.
#' @param mubeta k by m matrix. The first row represents the intercept (mu)
#'  for m components, and the second to k-th rows represent the
#'  coefficients on x for m components.
#' @param sigma m by 1 vector that represents sd of components.
#' @param x n by k-1 matrix that does NOT include a constant.
#' @return n by 1 vector that is formed by regressor x.
rnormregmix <- function (n, alpha, mubeta, sigma, x = NULL) {
  # Generates mixed normal random variables with regressor x
  # Input
  #  n : number of observations
  #  alpha  : m-vector
  #  mubeta  : k by m matrix
  #  sigma  : m-vector
  #  x : (n by k-1) matrix NOT including a constant
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



#' @description Computes omega_{j|i} defined in (2.1) of Maitra and Melnykov (2010).
#' @export
#' @title omega.ji
#' @name omega.ji
#' @param phi_i 3 by 1 column consisting of alpha, mu, sigma of ith component
#' @param phi_j 3 by 1 column consisting of alpha, mu, sigma of jth component
#' @return omega_{j|i}
#' @references Maitra, R., and Melnykov, V. (2010)
#' Simulating Data to Study Performance of Finite Mixture Modeling and Model-Based Clustering Algorithms,
#' \emph{Journal of Computational and Graphical Statistica},
#' \bold{19}, 354--376.
# Returns a misclassification rate omega_ji given two components i, j,
# i.e. the probability of choosing component j where
# the true model is ith component.
omega.ji <- function(phi_i, phi_j) {
  alpha_i <- phi_i[1]
  alpha_j <- phi_j[1]
  mu_i <- phi_i[2]
  mu_j <- phi_j[2]
  sigma_i <- phi_i[3]
  sigma_j <- phi_j[3]

  a <- (1/sigma_j^2 - 1/sigma_i^2)
  b <- mu_i / sigma_i^2 - mu_j / sigma_j^2
  c <- mu_j^2 / sigma_j^2 - mu_i^2 / sigma_i^2

  if (sigma_i == sigma_j)
    if (mu_i > mu_j)
      omega_ji = pnorm((2 * log(alpha_j/alpha_i) - c)/(2*b),
                       mean = mu_i, sd = sigma_i)
  else
    omega_ji = 1 - pnorm((2 * log(alpha_j/alpha_i) - c)/(2*b),
                         mean = mu_i, sd = sigma_i)
  else {
    d <- 2 * log(alpha_j * sigma_i / (alpha_i * sigma_j)) - c + (b^2 / a)
    da <- max(d/a, 0)
    if (sigma_i > sigma_j)
      omega_ji = pnorm(sqrt(da)-b/a, mean = mu_i, sd = sigma_i) -
        pnorm(-sqrt(da)-b/a, mean = mu_i, sd = sigma_i)
    else
      omega_ji = 1 +
      pnorm(-sqrt(da)-b/a, mean = mu_i, sd = sigma_i) -
      pnorm(sqrt(da)-b/a, mean = mu_i, sd = sigma_i)
  }

  omega_ji <- unname(omega_ji)
  return (omega_ji)
}

#' @description Computes \eqn{\omega_{12}}, the average misclassification rate between two normal mixtures defined in Maitra and Melnykov (2010).
#' @export
#' @title omega.12
#' @name omega.12
#' @param parlist parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gamma = (gamma_1, ..., gamma_m)).
#' @return misclassification rate \eqn{\omega_{12}}.
#' @references Maitra, R., and Melnykov, V. (2010)
#' Simulating Data to Study Performance of Finite Mixture Modeling and Model-Based Clustering Algorithms,
#' \emph{Journal of Computational and Graphical Statistica},
#' \bold{19}, 354--376.
omega.12 <- function(parlist)
  # Computes omega_{12} for testing H_0:m=2 against H_1:m=3
{
  phi1 <- c(alpha = parlist$alpha[1], mu = parlist$mu[1], sigma = parlist$sigma[1])
  phi2 <- c(alpha = parlist$alpha[2], mu = parlist$mu[2], sigma = parlist$sigma[2])

  part1 <- omega.ji(phi1, phi2)
  part2 <- omega.ji(phi2, phi1)
  
  omega.12 <- (part1 + part2) / 2

  return (omega.12)
}  # end function omega.12


#' Computes \eqn{\omega_{12}} and \eqn{\omega_{23}} defined in Maitra and Melnykov (2010).
#' @export
#' @title omega.123
#' @name omega.123
#' @param parlist parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gamma_1, ..., gamma_m)).
#' @return 2 by 1 vector whose first element is \eqn{\omega_{12}} and second element is \eqn{\omega_{23}}.
#' @references Maitra, R., and Melnykov, V. (2010)
#' Simulating Data to Study Performance of Finite Mixture Modeling and Model-Based Clustering Algorithms,
#' \emph{Journal of Computational and Graphical Statistica},
#' \bold{19}, 354--376.
omega.123 <- function(parlist)
{
  phi1 <- c(alpha = parlist$alpha[1], mu = parlist$mu[1], sigma = parlist$sigma[1])
  phi2 <- c(alpha = parlist$alpha[2], mu = parlist$mu[2], sigma = parlist$sigma[2])
  phi3 <- c(alpha = parlist$alpha[3], mu = parlist$mu[3], sigma = parlist$sigma[3])

  part1 <- omega.ji(phi1, phi2)
  part2 <- omega.ji(phi2, phi1)
  w12 <- (part1 + part2)/2

  part3 <- omega.ji(phi2, phi3)
  part4 <- omega.ji(phi3, phi2)
  w23 <- (part3 + part4)/2

  return(c(w12, w23))

}  # end function omega.123

#' @description Computes \eqn{\omega_{12}}, \eqn{\omega_{23}}, and \eqn{\omega_{34}} defined in Maitra and Melnykov (2010).
#' @export
#' @title omega.1234
#' @name omega.1234
#' @param parlist parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gamma_1, ..., gamma_m)).
#' @return 3 by 1 vector consisting of \eqn{\omega_{12}}, \eqn{\omega_{23}}, and \eqn{\omega_{34}}.
#' @references Maitra, R., and Melnykov, V. (2010)
#' Simulating Data to Study Performance of Finite Mixture Modeling and Model-Based Clustering Algorithms,
#' \emph{Journal of Computational and Graphical Statistica},
#' \bold{19}, 354--376.
omega.1234 <- function(parlist)
{
  phi1 <- c(alpha = parlist$alpha[1], mu = parlist$mu[1], sigma = parlist$sigma[1])
  phi2 <- c(alpha = parlist$alpha[2], mu = parlist$mu[2], sigma = parlist$sigma[2])
  phi3 <- c(alpha = parlist$alpha[3], mu = parlist$mu[3], sigma = parlist$sigma[3])
  phi4 <- c(alpha = parlist$alpha[4], mu = parlist$mu[4], sigma = parlist$sigma[4])

  part1 <- omega.ji(phi1, phi2)
  part2 <- omega.ji(phi2, phi1)
  w12 <- (part1 + part2)/2

  part3 <- omega.ji(phi2, phi3)
  part4 <- omega.ji(phi3, phi2)
  w23 <- (part3 + part4)/2

  part5 <- omega.ji(phi3, phi4)
  part6 <- omega.ji(phi4, phi3)
  w34 <- (part5 + part6)/2

  return(c(w12, w23, w34))

}  # end function omega.1234

