#' Computes the variance-covariance matrix of the MLE of m-component normal mixture.
#' @export
#' @title regmixVcov
#' @name regmixVcov
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x
#' @param coefficients (alpha_1, ..., alpha_m, mu_1, ..., mu_m, sigma_1, ..., sigma_m, gamma)
#' @param z n by p matrix of regressor associated with gamma
#' @param p Dimension of z
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
  #  coefficients : (alpha_1,...,alpha_m,mubeta_1^T, ...,mubeta_m^T,sigma_1, ..., sigma_m,gamma^T)
  #  z  : n by p matrix of regressor associated with gamma
  # Output
  #  vcov: variance-covariance matrix
  
  y     <- as.vector(y)
  n     <- length(y)
  len   <- length(coefficients)
  p     <- 0
  gamma <- NULL
  vcov.method <- match.arg(vcov.method)
  
  if (!is.null(z)) {
    z <- as.matrix(z)
    p <- ncol(z)
    gamma <- coefficients[(len-p+1):len]
  }
  
  x   <- as.matrix(x)
  x1  <- cbind(1, x)
  q1   <- ncol(x1)
  q   <- ncol(x)
  
  m  <- (len-p)/(3+q)
  if (round(m) != m) {
    stop("The dimension of the coefficients is incompatible with x and z. Please check the data.")
  }
  
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
    jpvt    <- integer(max(q1,p))    # pivots used in dgelsy
    
    if (is.null(z)) {
      setting <- c(n, m, q1, 1, 1, jpvt)
      out.p <- .C("regmixpmle", as.integer(setting), as.double(y), as.double(x1),
                  alphaset = as.double(alpha), mubetaset = as.double(mubeta), sigmaset = as.double(sigma),
                  as.double(sigma0), as.double(mu0), as.double(an), tau, as.integer(h), k, lub = double(2*m),
                  double(3*m), post = double(n*m),
                  loglikset = double(1), penloglikset = double(1),
                  notcg = integer(1), as.double(epsilon), double(n*(q1+1)), package = "normalregMix")
    } else {
      setting.z <- c(n, m, q1, p, 1, 1, jpvt)
      out.p <- .C("regmixpmle_z", as.integer(setting.z), as.double(y), as.double(x1), as.double(z),
                  alphaset = as.double(alpha), mubetaset = as.double(mubeta), sigmaset = as.double(sigma), gammaset = as.double(gamma),
                  as.double(sigma0), as.double(mu0), as.double(an), tau, as.integer(h), k, lub = double(2*m),
                  double(3*m), post = double(n*m),
                  loglikset = double(1), penloglikset = double(1),
                  notcg = integer(1), as.double(epsilon), double(n*(q1+3+p)), package = "normalregMix")
      # Adjust y
      y <- as.vector(y - z %*% gamma)
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
        Q.gamma.theta <- matrix(0, nrow=p, ncol=(q1+1)*m)  # p by (q1+1)*m matrix
        for (i in 1:m) {
          C.i <- array(C0[, i, 1, ], dim=c(n, q1+1))  # n by q1+1
          Q.i.1 <- colSums(tKR(-C.i+b[, i]*c0[, i, ], z*post[, i])) # p*(q1+1) vector
          Q.i.2 <- colSums(tKR(c0[, i, ]*post[, i], dbar))  # p*(q1+1) vector
          Q.gamma.theta[, ((q1+1)*(i-1)+1):((q1+1)*i)] <- matrix(Q.i.1+Q.i.2, nrow=p, ncol=q1+1)
        }
        
        Q.gamma.theta <- Q.gamma.theta[, c(p1, p2), drop=FALSE]  # p by (q1+1)*m
        w1 <- (post*b)%*%t(a) - rowSums(post*b)*t(abar)  # n by m-1
        Q.pi.gamma.0 <- colSums(tKR(w1, z))  # (m-1)*p vector
        Q.pi.gamma  <- matrix(Q.pi.gamma.0, nrow=m-1, ncol=p)
        Q.gamma     <- - t(z)%*%(z*rowSums(post*B)) -
          matrix(colSums(tKR(dbar, dbar)), nrow=p, ncol=p)
        
        I <- cbind(I, -rbind(Q.pi.gamma, t(Q.gamma.theta)))
        I <- rbind(I, -cbind(t(Q.pi.gamma), Q.gamma.theta, Q.gamma))
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

#' Computes the bootstrap critical values of the modified EM test.
#' @export
#' @title regmixCritBoot
#' @name regmixCritBoot
#' @param y n by 1 vector of data for y
#' @param x n by q vector of data for x 
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma 
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m), 
#' sigma = (sigma_1, ..., sigma_m), gamma = (gamma_1, ..., gamma_m))
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
  if (test.on) # initial values controlled by test.on
    set.seed(test.seed)
  
  y  <- as.vector(y)
  n  <- length(y)
  x  <- as.matrix(x)
  
  alpha   <- parlist$alpha
  mubeta  <- parlist$mubeta
  sigma   <- parlist$sigma
  gamma   <- parlist$gamma
  m       <- length(alpha)
  
  pvals <- NULL
  
  # Generate bootstrap observations
  ybset <- replicate(nbtsp, rnormregmix(n = n, x = x, alpha = alpha, mubeta = mubeta, sigma = sigma))
  
  if (!is.null(z)) {
    zgamma <- as.matrix(z) %*% gamma
    ybset <- ybset + as.vector(zgamma)
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


# regmixMEMtestNocrit <- function (y, x, m = 2, z = NULL, an = 2.2, ninits = 100) {
# # Compute the modified EM test statistic for testing H_0 of m components
# # against H_1 of m+1 components for a univariate finite mixture of normals
#
# regmix.pmle.result    <- regmixPMLE(y=y, x=x, m=m, z=z, vcov.method="none", ninits=ninits)
# loglik0 <- regmix.pmle.result$loglik
# par1    <- regmixMaxPhi(y=y, x=x, regmix.pmle.result=regmix.pmle.result$parlist, z=z, an=an, ninits=ninits)
# emstat  <- 2*(par1$loglik-loglik0)
#
# a <- list(emstat=emstat, parlist=regmix.pmle.result$parlist)
#
# a
#
# }  # end regmixMEMtestNocrit

#' Sequentially performs MEM test given the data for y and x on the null hypothesis H_0: m = m_0
#' where m_0 is in {1, 2, ..., maxm}
#' @export
#' @title regmixMEMtestSeq
#' @name regmixMEMtestSeq
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x
#' @param z n by p matrix of regressor associated with gamma
#' @param maxm The maximum number of components set as null hypothesis in the mixture
#' @param ninits The number of randomly drawn initial values.
#' @param maxit The maximum number of iterations.
#' @param nbtsp The number of bootstrap observations; by default, it is set to be 199
#' @param parallel Determines whether package \code{doParallel} is used for calculation
#' @param cl Cluster used for parallelization; if it is \code{NULL}, the system will automatically
#' create a new one for computation accordingly.
#' @return A list of with the following items:
#' \item{alpha}{maxm by maxm matrix, whose i-th column is a vector of alphas estimated given the null hypothesis m_0 = i}
#' \item{mu}{maxm by maxm matrix, whose i-th column is a vector of mus estimated given the null hypothesis m_0 = i}
#' \item{sigma}{maxm by maxm matrix, whose i-th column is a vector of sigmas estimated given the null hypothesis m_0 = i}
#' \item{beta}{A list of length maxm, whose i-th element is a q times i matrix of betas estimated given the null hypothesis m_0 = i}
#' \item{gamma}{maxm by maxm matrix, whose i-th column is a vector of gammas estimated given the null hypothesis m_0 = i}
#' \item{emstat}{A maxm vector of values of modified EM statistics of the model at m_0 = 1, 2, ..., maxm}
#' \item{pvals}{A maxm by 3 matrix whose i-th row indicates a vector of p-values at k = 1, 2, 3}
#' \item{aic}{A maxm vector of Akaike Information Criterion of the fitted model at m_0 = 1, 2, ..., maxm}
#' \item{bic}{A maxm vector of Bayesian Information Criterion of the fitted model at m_0 = 1, 2, ..., maxm}
#' \item{loglik}{A maxm vector of log-likelihood values of the model at m_0 = 1, 2, ..., maxm}
#' \item{penloglik}{A maxm vector of penalized log-likelihood values of the model at m_0 = 1, 2, ..., maxm}
#' @examples 
#' data(faithful)
#' attach(faithful)
#' regmixMEMtestSeq(y = eruptions, x = waiting)
regmixMEMtestSeq <- function (y, x, z = NULL, maxm = 3, ninits = 10, maxit = 2000,
                              nbtsp = 199, parallel = FALSE, cl = NULL) {
  # Compute the modified EM test statistic for testing H_0 of m components
  # against H_1 of m+1 components for a univariate finite mixture of normals
  
  y   <- as.vector(y)
  x		<- as.matrix(x)
  n   <- length(y)
  p   <- 0
  q   <- ncol(x)

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
  beta <- list()
  # Test H_0:m=1, H_0:m=2, ...
  
  for (m in 1:maxm){
    pmle.result   <- regmixPMLE(y = y, x = x, m = m, z = z, vcov.method = "none",
                                ninits = ninits, maxit = maxit)
    loglik[m] <- loglik0 <- pmle.result$loglik
    penloglik[m] <- penloglik0 <- pmle.result$penloglik
    aic[m]  <- pmle.result$aic
    bic[m]  <- pmle.result$bic
    
    parlist <- pmle.result$parlist
    alpha0  <- parlist$alpha
    mubeta0 <- parlist$mubeta
    mu0			<- mubeta0[1,]
    beta0   <- mubeta0[-1,]
    sigma0  <- parlist$sigma
    gamma0  <- parlist$gamma
    
    alpha[,m] <- c(alpha0, double(maxm - m))
    mu[,m]    <- c(mu0, double(maxm - m))
    sigma[,m] <- c(sigma0, double(maxm - m))
    beta[m]  <- c(beta0, double(maxm - m))
    
    cat(sprintf("%d-component model estimate:\n",m))
    
    
    tab = as.matrix(rbind(alpha0, mu0, sigma0, as.matrix(beta0)))
    
    print(tab)
    print(beta0)
    print(mubeta0)
    rownames(tab) <- c("alpha", "mu", "sigma", paste("beta", ".", 1:q, sep = ""))
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
      
      #     print(pmle.result$parlist)
      an    <- anFormula(parlist = parlist, m = m, n = n)
      #     print(an)
      par1  <- regmixMaxPhi(y = y, x = x, parlist = parlist, z = z, an = an, 
                            ninits = ninits, maxit = maxit, parallel = parallel)
      #     print(loglik0)
      #     print(par1)
      emstat.m  <- 2*(par1$penloglik - loglik0)
      
      cat(c("modified EM-test statitic ", sprintf('%.3f',emstat.m)),"\n")
      if (m <=3 ) {
        em.out <- regmixCrit(y = y, x = x, parlist = parlist, z = z, 
                             values = emstat.m)
        cat(c("asymptotic p-values       ", sprintf('%.3f',em.out$pvals)),"\n \n")
      } else {
        em.out <- regmixCritBoot(y = y, x = x, parlist=parlist, z = z, 
                                 values = emstat.m, ninits = ninits, nbtsp = nbtsp, 
                                 parallel = parallel, cl = cl)
        cat(c("bootstrap p-values        ", sprintf('%.3f',em.out$pvals)),"\n \n")
      }
      # noncg.rate[m]   <- par1$noncg.rate
      pvals[m,]     <- em.out$pvals
      emstat[m,]    <- emstat.m
    }
  }
  a = list(alpha = alpha, mu = mu, sigma = sigma, beta = beta, gamma = gamma, 
           emstat = emstat, pvals = pvals, aic = aic, bic = bic, 
           loglik = loglik, penloglik = penloglik)
  
  a
  
}  # end regmixMEMtestSeq

#' Performs MEM test given the data for y and x on the null hypothesis H_0: m = m_0.
#' @export
#' @title regmixMEMtest
#' @name regmixMEMtest
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x
#' @param m The number of components in the mixture defined by a null hypothesis, m_0
#' @param z n by p matrix of regressor associated with gamma
#' @param tauset A set of initial tau value candidates
#' @param an a term used for penalty function
#' @param ninits The number of randomly drawn initial values.
#' @param crit.method Method used to compute the variance-covariance matrix, one of \code{"none"}, 
#' \code{"asy"}, and \code{"boot"}. The default option is \code{"none"}. When \code{method = "asy"}, 
#' the p-values are computed based on an asymptotic method. When \code{method = "OPG"}, 
#' the p-values are generated by bootstrapping.
#' @param nbtsp The number of bootstrap observations; by default, it is set to be 199
#' @param cl Cluster used for parallelization; if it is \code{NULL}, the system will automatically
#' create a new one for computation accordingly.
#' @param parallel Determines whether package \code{doParallel} is used for calculation
#' @return A list of class \code{normalMix} with items:
#' \item{coefficients}{A vector of parameter estimates. Ordered as \eqn{\alpha_1,\ldots,\alpha_m,\mu_1,\ldots,\mu_m,\sigma_1,\ldots,\sigma_m,\gamma}.}
#' \item{parlist}{The parameter estimates as a list containing alpha, mu, and sigma (and gamma if z is included in the model).}
#' \item{vcov}{The estimated variance-covariance matrix.}
#' \item{loglik}{The maximized value of the log-likelihood.}
#' \item{penloglik}{The maximized value of the penalized log-likelihood.}
#' \item{aic}{Akaike Information Criterion of the fitted model.}
#' \item{bic}{Bayesian Information Criterion of the fitted model.}
#' \item{postprobs}{n by m matrix of posterior probabilities for observations}
#' \item{indices}{n by 1 vector of integers that indicates the indices of components 
#' each observation belongs to based on computed posterior probabilities}
#' \item{call}{The matched call.}
#' \item{m}{The number of components in the mixture.}
#' @examples 
#' data(faithful)
#' attach(faithful)
#' regmixMEMtest(y = eruptions, x = waiting, m = 1, crit.method = "asy")
#' regmixMEMtest(y = eruptions, x = waiting, m = 2)
regmixMEMtest <- function (y, x, m = 2, z = NULL, tauset = c(0.1,0.3,0.5), 
                           an = NULL, ninits = 100,
                           crit.method = c("none", "asy", "boot"), nbtsp = 199,
                           cl = NULL, 
                           parallel = FALSE) {
  # Compute the modified EM test statistic for testing H_0 of m components
  # against H_1 of m+1 components for a univariate finite mixture of normals
  y <- as.vector(y)
  x <- as.matrix(x)
  q <- ncol(x)
  
  if (!is.null(z)) 
    z <- as.matrix(z)
  
  regmix.pmle.result    <- regmixPMLE(y=y, x=x, m=m, z=z, vcov.method="none", ninits=ninits)
  loglik0 <- regmix.pmle.result$loglik
  
  if (is.null(an)) 
    an <- anFormula(parlist = regmix.pmle.result$parlist, m = m, n = n, q = q)
  
  
  par1    <- regmixMaxPhi(y=y, x=x, parlist=regmix.pmle.result$parlist, z=z, 
                          an=an, tauset = tauset, ninits=ninits,
                          parallel = parallel, cl = cl)
  # use the penalized log-likelihood.
  emstat  <- 2*(par1$penloglik-loglik0)
  
  if (crit.method == "asy"){
    result  <- regmixCrit(y=y, x=x, parlist=regmix.pmle.result$parlist, z=z, values=emstat,
                          parallel=parallel, cl=cl, nrep=1000, ninits.crit=25)
  } else if (crit.method == "boot") {
    result  <- regmixCritBoot(y=y, x=x, parlist=regmix.pmle.result$parlist, z=z, values=emstat,
                              ninits=ninits, nbtsp=nbtsp, parallel=parallel, cl=cl)
  } else {
    result <- list()
    result$crit <- result$pvals <- rep(NA,3)
  }
  
  a <- list(emstat = emstat, pvals = result$pvals, crit = result$crit, parlist = regmix.pmle.result$parlist,
            call = match.call(), m = m, crit.method = crit.method, nbtsp = nbtsp,
            label = "MEMtest")
  
  class(a) <- "normalregMix"
  
  a
  
}  # end regmixMEMtest

## TODO: verb is not used.
#' Compute ordinary & penalized log-likelihood ratio resulting from MEM algorithm at k=1,2,3.
#' @export
#' @title regmixMaxPhi
#' @name regmixMaxPhi
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma 
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m), 
#' sigma = (sigma_1, ..., sigma_m), gamma = (gamma_1, ..., gamma_m))
#' @param z n by p matrix of regressor associated with gamma
#' @param an a term used for penalty function
#' @param tauset A set of initial tau value candidates
#' @param ninits The number of randomly drawn initial values.
#' @param epsilon.short The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param epsilon The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}. 
#' @param maxit.short The maximum number of iterations in short EM.
#' @param maxit The maximum number of iterations.
#' @param parallel Determines whether package \code{doParallel} is used for calculation
#' @param cl Cluster used for parallelization; if it is \code{NULL}, the system will automatically
#' create a new one for computation accordingly. 
#' @return A list with items:
#' \item{loglik}{Log-likelihood resulting from MEM algorithm at k=1,2,3.}
#' \item{penloglik}{Penalized log-likelihood resulting from MEM algorithm at k=1,2,3.}
regmixMaxPhi <- function (y, x, parlist, z = NULL, an, tauset = c(0.1,0.3,0.5), 
                          ninits = 100,
                          epsilon.short = 1e-02, epsilon = 1e-08,
                          maxit.short = 500, maxit = 2000,
                          verb = FALSE, 
                          parallel = FALSE,
                          cl = NULL) {
  # Given a parameter estiamte of an m component model and tuning paramter an,
  # maximize the objective function for computing the modified EM test statistic
  # for testing H_0 of m components against H_1 of m+1 for a univariate finite mixture of normals
  
  warn  <- options(warn=-1) # Turn off warnings

  # q1 = dim(X) + dim(mu) = q + 1
  q1 <- ncol(x) + 1
  p <- 0
  m <- length(parlist$alpha)
  
  if (!is.null(z)) 
    p <- ncol(z)
  jpvt  <- integer(max(q1,p))    # pivots used in dgelsy
  
  ninits.short <- ninits*4*(q1+p)*m
  
  loglik.all <- matrix(0,nrow=m*length(tauset),ncol=3)
  penloglik.all <- matrix(0,nrow=m*length(tauset),ncol=3)
  
  if (parallel) {
    if (is.null(cl))
      cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    results <- foreach (t = 1:length(tauset),
                        .export = 'regmixPhiStep', .combine = c)  %:%
      foreach (h = 1:m) %dopar% {
        regmixPhiStep (c(h, tauset[t]), y, x, parlist, z = NULL, p, jpvt,
                       an,
                       ninits, ninits.short,
                       epsilon.short, epsilon,
                       maxit.short, maxit,
                       verb) }
    on.exit(cl)
    loglik.all <- t(sapply(results, "[[", "loglik"))
    penloglik.all <- t(sapply(results, "[[", "penloglik"))
  }
  else
    for (h in 1:m) 
      for (t in 1:length(tauset)) {
        rowindex <- (t-1)*m + h
        tau <- tauset[t]
        result <- regmixPhiStep(c(h, tau), y, x, parlist, z = NULL, p, jpvt, 
                                an,
                                ninits, ninits.short,
                                epsilon.short, epsilon,
                                maxit.short, maxit,
                                verb)
        loglik.all[rowindex,] <- result$loglik
        penloglik.all[rowindex,] <- result$penloglik
      }
  
  loglik <- apply(loglik.all, 2, max)  # 3 by 1 vector
  penloglik <- apply(penloglik.all, 2, max)  # 3 by 1 vector
  
  out <- list(loglik = loglik, penloglik = penloglik)
  
  out
  
}  # end regmixMaxPhi

#' Given a pair of h and tau and data, compute ordinary & 
#' penalized log-likelihood ratio resulting from MEM algorithm at k=1,2,3, tailored for parallelization.
#' @export
#' @title regmixPhiStep
#' @name regmixPhiStep
#' @param htaupair A set of h and tau
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma 
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m), 
#' sigma = (sigma_1, ..., sigma_m), gamma = (gamma_1, ..., gamma_m))
#' @param z n by p matrix of regressor associated with gamma
#' @param p Dimension of z
#' @param jpvt Pivots used in dgelsy
#' @param an a term used for penalty function
#' @param ninits The number of randomly drawn initial values.
#' @param ninits.short The number of candidates used to generate an initial phi, in short MEM
#' @param epsilon.short The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param epsilon The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}. 
#' @param maxit.short The maximum number of iterations in short EM.
#' @param maxit The maximum number of iterations.
#' @return A list of phi, log-likelihood, and penalized log-likelihood resulting from MEM algorithm.
regmixPhiStep<- function (htaupair, y, x, parlist, z = NULL, p, jpvt, 
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
  mu0h <- c(0,mubeta0[1,],0)        # m+1 by 1
  sigma0 <- parlist$sigma
  sigma0h <- c(sigma0[1:h],sigma0[h:m])  # m+1 by 1
  gamma0 <- parlist$gamma
  
  # generate initial values
  tmp <- regmixPhiInit(y = y, x = x, z = z, parlist=parlist, h=h, tau, ninits = ninits.short)
  # tmp$alpha, tmp$sigma: m+1 by ninits matrix, tmp$mubeta: p by m+1 by ninits array
  # tmp$gamma: p by ninits matrix
  
  alphaset.s   <- tmp$alpha
  mubetaset.s <- matrix(tmp$mubeta, nrow=m1*q1)
  sigmaset.s   <- tmp$sigma
  gammaset.s  <- tmp$gamma
  
  if (is.null(z)) {
    setting <- c(n,m1,q1,ninits.short,maxit.short,jpvt)  # configulation parameters
    out.short <- .C("regmixpmle", as.integer(setting), as.double(y), as.double(x1),
                    alphaset = as.double(alphaset.s), mubetaset = as.double(mubetaset.s), sigmaset = as.double(sigmaset.s),
                    as.double(sigma0h), as.double(mu0h), as.double(an), 
                    as.double(tau), as.integer(h), as.integer(k), lub = double(2*m1),
                    double(3*m1), post = double(n*m1),
                    loglikset = double(ninits.short), penloglikset = double(ninits.short),
                    notcg = integer(ninits.short), as.double(epsilon.short), double(n*(q1+1)))
  } else {
    setting.z <- c(n,m1,q1,p,ninits.short,maxit.short,jpvt)  # configulation parameters
    out.short <- .C("regmixpmle_z", as.integer(setting.z), as.double(y), as.double(x1), as.double(z),
                    alphaset = as.double(alphaset.s), mubetaset = as.double(mubetaset.s), sigmaset = as.double(sigmaset.s), gammaset = as.double(gammaset.s),
                    as.double(sigma0h), as.double(mu0h), as.double(an), 
                    as.double(tau), as.integer(h), as.integer(k), lub = double(2*m1),
                    double(3*m1), post = double(n*m1),
                    loglikset = double(ninits.short), penloglikset = double(ninits.short),
                    notcg = integer(ninits.short), as.double(epsilon.short), double(n*(q1+3+p)))
  }
  
  # if (verb && any(out.short$notcg)) {
  #     cat(sprintf("non-convergence rate at short-EM = %.3f\n",mean(out.short$notcg)))
  # }
  
  penloglik.short <- out.short$penloglikset    # extract the 4th argument = penloglik
  oo <- order(penloglik.short, decreasing = TRUE)
  oo.inits <- oo[1:ninits]
  
  # long EM
  
  alphaset  <- alphaset.s[,oo.inits]
  mubetaset <- mubetaset.s[,oo.inits]
  sigmaset  <- sigmaset.s[,oo.inits]
  
  if (is.null(z)) {
    setting <- c(n,m1,q1,ninits,maxit,jpvt)
    out <- .C("regmixpmle", as.integer(setting), as.double(y), as.double(x1),
              alphaset = as.double(alphaset), mubetaset = as.double(mubetaset), sigmaset = as.double(sigmaset),
              as.double(sigma0h), as.double(mu0h), as.double(an), 
              as.double(tau), as.integer(h), as.integer(k), lub = double(2*m1),
              double(3*m1), post = double(n*m1),
              loglikset = double(ninits), penloglikset = double(ninits),
              notcg = integer(ninits), as.double(epsilon), double(n*(q1+1)))
  } else {
    gammaset  <- gammaset.s[,oo.inits]
    setting.z <- c(n,m1,q1,p,ninits,maxit,jpvt)
    out <- .C("regmixpmle_z", as.integer(setting.z), as.double(y), as.double(x1), as.double(z),
              alphaset = as.double(alphaset), mubetaset = as.double(mubetaset), sigmaset = as.double(sigmaset), gammaset = as.double(gammaset),
              as.double(sigma0h), as.double(mu0h), as.double(an), 
              as.double(tau), as.integer(h), as.integer(k), lub = double(2*m1),
              double(3*m1), post = double(n*m1),
              loglikset = double(ninits), penloglikset = double(ninits),
              notcg = integer(ninits), as.double(epsilon), double(n*(q1+3+p)))
  }
  
  if (any(out$notcg)) {
    cat(sprintf("non-convergence rate = %.3f\n",mean(out$notcg)))
  }
  
  index    <- which.max(out$penloglikset)
  alpha   <- out$alphaset[(m1*(index-1)+1):(m1*index)]
  mubeta  <- matrix(out$mubetaset[(q1*m1*(index-1)+1):(q1*m1*index)], nrow=q1)
  sigma   <- out$sigmaset[(m1*(index-1)+1):(m1*index)]
  gamma    <- out$gammaset[(p*(index-1)+1):(p*index)]
  penloglik <- out$penloglikset[index]
  loglik     <- out$loglikset[index]
  
  mu.order  <- order(mubeta[1,])
  alpha     <- alpha[mu.order]
  mubeta    <- mubeta[,mu.order]
  sigma      <- sigma[mu.order]
  
  a0 = list(alpha = alpha, mubeta = mubeta, sigma = sigma, gamma = gamma,
            penloglik = penloglik, loglik=loglik,
            ft="normalmixEM")
  
  a     <- regmixPhi2(y, x, z, a0, sigma0 = sigma0, h=h, tau=tau, an=an)
  a[[1]]  <- list(alpha = alpha, mubeta = mubeta, sigma = sigma, gamma = gamma,
                  penloglik = a0$penloglik, loglik=a0$loglik)
  
  phi1 <- c(a[[1]], h = h)
  loglik <- sapply(a, "[[", "loglik")    # extract loglik at q1=1,2,3
  penloglik <- sapply(a, "[[", "penloglik")    # extract penloglik at q1=1,2,3
  
  return (list(phi1 = phi1, loglik = loglik, penloglik = penloglik))
  
}

#' Starting from a parameter estimate, take two EM steps to compute the
#' local modified EM test statistic for testing H_0 of m components
#' against H_1 of m+1 at K=2,3 for a univariate normal finite mixture.
#' @export
#' @title regmixPhi2
#' @name regmixPhi2
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x
#' @param z n by p matrix of regressor associated with gamma
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma 
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m), 
#' sigma = (sigma_1, ..., sigma_m), gamma = (gamma_1, ..., gamma_m))
#' @param sigma0 m by 1 vector of sigma 
#' @param h h used as index for pivoting
#' @param tau Tau used to split the h-th component
#' @param an a term used for penalty function
#' @return A list of length 2, where the first element is the local modified 
#' MEM test statistic at k=2, the second element is 
#' the local modified MEM test statistic at k=3
regmixPhi2 <- function (y, x, z=NULL, parlist, sigma0, h, tau, an) {
  # Starting from a parameter estimate, take two EM steps to compute the
  # local modified EM test statistic for testing H_0 of m components
  # against H_1 of m+1 at q1=2,3 for a finite mixture of regressions
  # In input,
  #  y    : n by 1 vector of dependent variable
  #  x    : n by q1-1 matrix of regressor NOT including an intercept
  #  z    : n by p matrix of regressor associated with gamma
  
  n     <- length(y)
  x1     <- cbind(1,x)
  q1     <- ncol(x1)
  alpha  <- parlist$alpha
  mubeta  <- parlist$mubeta
  sigma  <- parlist$sigma
  gamma  <- parlist$gamma
  m1     <- length(alpha) # m1 = m+1
  m     <- m1-1
  sigma0h <- c(sigma0[1:h],sigma0[h:m])
  a     <- vector('list',3)
  mu0   <- double(m1)
  # In the old code, it used h0 as a parameter for the C code instead of h 
  # which is actually taken as a parameter in this function 
  # (in fact, h is used nowhere except configuring sigma0h)
  # h0     <- 0
  
  p      <- 0
  if (!is.null(z)) {p <- ncol(z)}
  
  jpvt    <- integer(max(q1,p))
  setting <- c(n,m1,q1,1,1,jpvt)
  setting.z <- c(n,m1,q1,p,1,1,jpvt)
  
  for (k in 2:3) {
    a[[k]] <- list(penloglik = -Inf, loglik = -Inf)
  }
  
  for (k in 2:3) {
    # Two EM steps
    if (is.null(z)) {
      out <- .C("regmixpmle", as.integer(setting), as.double(y), as.double(x1),
                alpha = as.double(alpha), mubeta = as.double(mubeta), 
                sigma = as.double(sigma),
                as.double(sigma0h), as.double(mu0), as.double(an), 
                as.double(tau), as.integer(h), as.integer(k), lub = double(2*m1),
                double(3*m1), post = double(n*m1),
                loglik = double(1), penloglik = double(1),
                notcg = integer(1), tol = as.double(1e-8), double(n*(q1+1)))
    } else {
      out <- .C("regmixpmle_z", as.integer(setting.z), as.double(y), as.double(x1), as.double(z),
                alpha = as.double(alpha), mubeta = as.double(mubeta), 
                sigma = as.double(sigma), gamma = as.double(gamma),
                as.double(sigma0h), as.double(mu0), as.double(an), 
                as.double(tau), as.integer(h), as.integer(k), lub = double(2*m1),
                double(3*m1), post = double(n*m1),
                loglik = double(1), penloglik = double(1),
                notcg = integer(1), tol = as.double(1e-8), double(n*(q1+3+p)))
    }
    
    # Check singularity: if singular, break from the loop
    if ( any(out$sigma < 1e-06) || any(out$alpha < 1e-06) || is.na(sum(out$alpha)) ) {
      break
    }
    
    alpha    <- out$alpha
    mubeta  <- out$mubeta
    sigma    <- out$sigma
    gamma    <- out$gamma
    loglik    <- out$loglik
    penloglik <- out$penloglik
    
    a[[k]] <- list(alpha = alpha, mubeta = mubeta, sigma = sigma, gamma = gamma, penloglik = penloglik, loglik = loglik)
    
  }
  
  a
  
}  # end function regmixPhi2

#' Generates lists of parameters for initial candidates used by 
#' the modified EM test for mixture of normals.
#' @export
#' @title regmixPhiInit
#' @name regmixPhiInit
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x 
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma 
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m), 
#' sigma = (sigma_1, ..., sigma_m), gamma = (gamma_1, ..., gamma_m))
#' @param z n by p matrix of regressor associated with gamma
#' @param h h used as index for pivoting
#' @param tau Tau used to split the h-th component
#' @param ninits number of initial values to be generated
#' @return A list with the following items:
#' \item{alpha}{m+1 by ninits matrix for alpha}
#' \item{mubeta}{q+1 by m+1 by ninits array for mu and beta}
#' \item{sigma}{m+1 by ninits matrix for sigma}
#' \item{gamma}{m+1 by ninits matrix for gamma}
regmixPhiInit <- function (y, x, z = NULL, parlist, h, tau, ninits = 1)
{
  # Generate initial values used by the modified EM test for mixture of regressions
  # input is
  #  y    : n by 1 vector of dependent variable
  #  x    : n by q1-1 matrix of regressor NOT including an intercept
  #  z    : n by p matrix of regressor associated with gamma
  # parlist  : list of alpha, mubeta, sigma, gamma from an m-component model
  #  h    : h in the modified EM test
  #  ninits  : number of initial values to be generated
  # output is a list that contains
  #  alpha  : m+1 by ninits matrix
  #  mubeta  : q1 by m+1 by ninits array
  #  sigma  : m+1 by ninits matrix
  #  gamma  : p by ninits matrix
  if (test.on) # initial values controlled by test.on
    set.seed(test.seed)
  
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
  gamma  <- NULL
  
  if (!is.null(z)) {
    gamma0  <- parlist$gamma
    p     <- ncol(z)
    y     <- as.vector(y - z %*% gamma0)
    gamma1   <- c(rep(1,p),runif(p*(ninits-1),min=-2,max=2))*gamma0
    gamma   <- matrix(gamma1,nrow=p)
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
  
  mu <- c(mu.hyp,runif((m+1)*(ninits-1),min=lb,max=ub))
  beta1 <- c(rep(1,(m+1)*(q1-1)),runif((m+1)*(q1-1)*(ninits-1),min=-2,max=2))*as.vector(beta.hyp)
  beta <- matrix(beta1,nrow=q1-1)
  
  mubeta <- rbind(mu,beta)
  mubeta <- array(mubeta,dim=c(q1,m+1,ninits))
  
  sigma.hyp <- c(sigma0[1:h],sigma0[h:m])  # m+1 by 1
  sigma1 <- c(rep(1,m+1),runif((m+1)*(ninits-1),min=0.25,max=2))*sigma.hyp
  sigma <- matrix(sigma1,nrow=m+1)
  
  alpha.hyp <- c(alpha0[1:h],alpha0[h:m])  # m+1 by 1
  alpha.hyp[h:(h+1)] <- c(alpha.hyp[h]*tau,alpha.hyp[h+1]*(1-tau))
  alpha <- matrix(rep.int(alpha.hyp,ninits),nrow=m+1)
  
  list(alpha = alpha, mubeta = mubeta, sigma = sigma, gamma = gamma)
  
}  # end function regmixPhiInit

#' Estimates parameters of a finite mixture of univariate normals by the method of penalized maximum likelhood.
#' @export
#' @title regmixPMLE
#' @name regmixPMLE
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x
#' @param m The number of components in the mixture
#' @param z n by p matrix of regressor associated with gamma
#' @param vcov.method Method used to compute the variance-covariance matrix, one of \code{"Hessian"} and \code{"OPG"}.
#' The default option is \code{"Hessian"}. When \code{method = "Hessian"}, the variance-covarince matrix is 
#' estimated by the Hessian using the formula given in Boldea and Magnus (2009). 
#' When \code{method = "OPG"}, the outer product of gradients is used.
#' @param ninits The number of randomly drawn initial values.
#' @param epsilon The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}. 
#' @param maxit The maximum number of iterations.
#' @param epsilon.short The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param maxit.short The maximum number of iterations in short EM.
#' @return  A list of class \code{normalMix} with items:
#' \item{coefficients}{A vector of parameter estimates. Ordered as \eqn{\alpha_1,\ldots,\alpha_m,\mu_1,\ldots,\mu_m,\sigma_1,\ldots,\sigma_m,\gamma}.}
#' \item{parlist}{The parameter estimates as a list containing alpha, mu, and sigma (and gamma if z is included in the model).}
#' \item{vcov}{The estimated variance-covariance matrix.}
#' \item{loglik}{The maximized value of the log-likelihood.}
#' \item{penloglik}{The maximized value of the penalized log-likelihood.}
#' \item{aic}{Akaike Information Criterion of the fitted model.}
#' \item{bic}{Bayesian Information Criterion of the fitted model.}
#' \item{postprobs}{n by m matrix of posterior probabilities for observations}
#' \item{indices}{n by 1 vector of integers that indicates the indices of components 
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
                        ninits = 100, epsilon = 1e-08, maxit = 2000,
                        epsilon.short = 1e-02, maxit.short = 500) {
  y   <- as.vector(y)
  x   <- as.matrix(x)   # n by (q1-1) matrix
  n   <- length(y)
  if (nrow(x) != n) { stop("y and x must have the same number of rows.") }
  x1  <- cbind(1, x)
  q1   <- ncol(x1)
  
  p       <- 0
  gamma   <- NULL
  ninits.short <- ninits*4*(q1+p)*m
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
    if (!is.null(z)) {gamma <- unname(ls.out$coeff[(q1+1):(q1+p)])}
    res     <- ls.out$residuals
    sigma   <- sqrt(mean(res*res))
    loglik  <- - (n/2)*(1 + log(2*pi) + 2*log(sigma))
    aic     <- -2*loglik + 2*npar
    bic     <- -2*loglik + log(n)*npar
    penloglik <- loglik
    
    parlist <- list(alpha = 1, mubeta = mubeta, sigma = sigma, gamma = gamma)
    coefficients <- c(alpha = 1, mubeta = mubeta, sigma = sigma, gamma = gamma)
    postprobs <- rep(1, n)
    
  } else {  # m >= 2
    
    # generate initial values
    tmp <- regmixPMLEinit(y = y, x = x, z = z, ninits = ninits.short, m = m)
    alphaset.s  <- tmp$alpha
    mubetaset.s <- matrix(tmp$mubeta, nrow=m*q1)
    sigmaset.s  <- tmp$sigma
    gammaset.s  <- tmp$gamma
    
    h       <- 0    # setting h=0 gives PMLE
    tau     <- 0.5  # setting tau=0.5 gives PMLE
    k <- 0 # setting k=0 gives PMLE
    
    sigma0  <- rep(sd0, m)
    mu0     <- double(m)    # dummy
    an      <- 1/n  # penalty term for variance
    
    jpvt <- integer(max(q1,p))  # pivots used in dgelsy
    
    # short EM
    if (is.null(z)) {
      setting <- c(n,m,q1,ninits.short,maxit.short,jpvt)  # configulation parameters
      out.short <- .C("regmixpmle", as.integer(setting), as.double(y), as.double(x1),
                      alphaset = as.double(alphaset.s), mubetaset = as.double(mubetaset.s), sigmaset = as.double(sigmaset.s),
                      as.double(sigma0), as.double(mu0), as.double(an), 
                      as.double(tau), as.integer(h), as.integer(k), 
                      lub = double(2*m),
                      double(3*m), post = double(n*m),
                      loglikset = double(ninits.short), penloglikset = double(ninits.short),
                      notcg = integer(ninits.short), as.double(epsilon.short), double(n*(q1+1)), package = "normalregMix")
    } else {
      setting.z <- c(n,m,q1,p,ninits.short,maxit.short,jpvt)  # configulation parameters
      out.short <- .C("regmixpmle_z", as.integer(setting.z), as.double(y), as.double(x1), as.double(z),
                      alphaset = as.double(alphaset.s), mubetaset = as.double(mubetaset.s), sigmaset = as.double(sigmaset.s), gammaset = as.double(gammaset.s),
                      as.double(sigma0), as.double(mu0), as.double(an), 
                      as.double(tau), as.integer(h), as.integer(k), 
                      lub = double(2*m),
                      double(3*m), post = double(n*m),
                      loglikset = double(ninits.short), penloglikset = double(ninits.short),
                      notcg = integer(ninits.short), as.double(epsilon.short), double(n*(q1+3+p)), package = "normalregMix")
    }
    
    penloglik.short <- out.short$penloglikset    # extract the 4th argument = penloglik
    oo <- order(penloglik.short, decreasing = TRUE)
    oo.inits <- oo[1:ninits]
    
    # long EM
    
    alphaset  <- alphaset.s[,oo.inits]
    mubetaset <- mubetaset.s[,oo.inits]
    sigmaset  <- sigmaset.s[,oo.inits]
    
    if (is.null(z)) {
      setting <- c(n,m,q1,ninits,maxit,jpvt)
      out <- .C("regmixpmle", as.integer(setting), as.double(y), as.double(x1),
                alphaset = as.double(alphaset), mubetaset = as.double(mubetaset), sigmaset = as.double(sigmaset),
                as.double(sigma0), as.double(mu0), as.double(an), 
                as.double(tau), as.integer(h), as.integer(k), 
                lub = double(2*m),
                double(3*m), post = double(n*m),
                loglikset = double(ninits), penloglikset = double(ninits),
                notcg = integer(ninits), as.double(epsilon), double(n*(q1+1)), package = "normalregMix")
    } else {
      setting.z <- c(n,m,q1,p,ninits,maxit,jpvt)
      gammaset  <- gammaset.s[,oo.inits]
      out <- .C("regmixpmle_z", as.integer(setting.z), as.double(y), as.double(x1), as.double(z),
                alphaset = as.double(alphaset), mubetaset = as.double(mubetaset), sigmaset = as.double(sigmaset), gammaset = as.double(gammaset),
                as.double(sigma0), as.double(mu0), as.double(an), 
                as.double(tau), as.integer(h), as.integer(k), 
                lub = double(2*m),
                double(3*m), post = double(n*m),
                loglikset = double(ninits), penloglikset = double(ninits),
                notcg = integer(ninits), as.double(epsilon), double(n*(q1+3+p)), package = "normalregMix")
    }
    
    #  if (mean(out$notcg) >= 0.9) {
    #  warning(sprintf("The EM algorithm failed to converge in %d%% of the initial values.
    # Try increasing maxit.", 100*mean(out$notcg)))
    #  }
    
    index     <- which.max(out$penloglikset)
    alpha     <- out$alphaset[(m*(index-1)+1):(m*index)]
    mubeta    <- matrix(out$mubetaset[(q1*m*(index-1)+1):(q1*m*index)], nrow=q1)
    sigma     <- out$sigmaset[(m*(index-1)+1):(m*index)]
    gamma     <- out$gammaset[(p*(index-1)+1):(p*index)]
    penloglik <- out$penloglikset[index]
    loglik    <- out$loglikset[index]
    
    # Compute posterior probabilities
    if (is.null(z)) {
      setting <- c(n,m,q1,1,1,jpvt)
      out.p <- .C("regmixpmle", as.integer(setting), as.double(y), as.double(x1),
                  alphaset = as.double(alpha), mubetaset = as.double(mubeta), 
                  sigmaset = as.double(sigma),
                  as.double(sigma0), as.double(mu0), as.double(an), 
                  as.double(tau), as.integer(h), as.integer(k), 
                  lub = double(2*m),
                  double(3*m), post = double(n*m),
                  loglikset = double(1), penloglikset = double(1),
                  notcg = integer(1), as.double(epsilon), double(n*(q1+1)), package = "normalregMix")
    } else {
      setting.z <- c(n,m,q1,p,1,1,jpvt)
      out.p <- .C("regmixpmle_z", as.integer(setting.z), as.double(y), as.double(x1), as.double(z),
                  alphaset = as.double(alpha), mubetaset = as.double(mubeta), 
                  sigmaset = as.double(sigma), gammaset = as.double(gamma),
                  as.double(sigma0), as.double(mu0), as.double(an), 
                  as.double(tau), as.integer(h), as.integer(k), 
                  lub = double(2*m),
                  double(3*m), post = double(n*m),
                  loglikset = double(1), penloglikset = double(1),
                  notcg = integer(1), as.double(epsilon), double(n*(q1+3+p)), package = "normalregMix")
    }
    
    aic <- -2*loglik + 2*npar
    bic <- -2*loglik + log(n)*npar
    
    mu.order  <- order(mubeta[1,])
    alpha     <- alpha[mu.order]
    mubeta    <- mubeta[,mu.order]
    sigma     <- sigma[mu.order]
    
    mu    <- mubeta[1,]
    beta  <- mubeta[-1,]
    
    postprobs <- matrix(out.p$post, nrow=n)
    postprobs <- postprobs[, mu.order]
    colnames(postprobs) <- c(paste("comp", ".", 1:m, sep = ""))
    
    mubeta.name <- matrix(0,nrow = q1, ncol = m)
    mubeta.name[1,] <- paste("mu", 1:m, sep = "")
    
    if (q1 == 2) {
      mubeta.name[2,] <- paste("beta", 1:m,  sep = "")
    } else {
      for (i in 1:(q1-1)) {
        for (j in 1:m) {
          mubeta.name[i+1,j] <- paste("beta", j, i, sep = "")
        }
      }    
    }
    
    parlist <- list(alpha = alpha, mubeta = mubeta, sigma = sigma, gamma = gamma)
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
            indices = apply(postprobs, 1, function(i) (which(i==max(i)))),
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
#' @param z n by p matrix of regressor associated with gamma
#' @param ninits number of initial values to be generated
#' @param m The number of components in the mixture
#' @return A list with the following items:
#' \item{alpha}{m by ninits matrix for alpha}
#' \item{mubeta}{q+1 by m by ninits array for mu and beta}
#' \item{sigma}{m by ninits matrix for sigma}
#' \item{gamma}{m by ninits matrix for gamma}
regmixPMLEinit <- function (y, x, z = NULL, ninits = 1, m = 2)
{  
  if (test.on) # initial values controlled by test.on
    set.seed(test.seed)
  
  n  <- length(y)
  q1  <- ncol(x)+1
  p  <- ncol(z)
  
  gamma <- NULL
  if (!is.null(z)) {
    out     <- lsfit(cbind(x, z), y)
    gamma0  <- out$coef[(q1+1):(q1+p)]
    gamma   <- matrix(runif(p*ninits, min=0.5, max=1.5), nrow=p)*gamma0
    mubeta_hat <- out$coef[1:q1]
    y     <- y - z %*% gamma0
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
  
  mu0   <- mubeta_hat[1]
  beta0 <- mubeta_hat[-1]
  minMU <- min(y - x %*% beta0)
  maxMU <- max(y - x %*% beta0)
  
  mu     <- runif(m*ninits, min=minMU, max=maxMU)
  beta   <- matrix(runif(m*(q1-1)*ninits, min=-2, max=2), nrow=q1-1)*beta0
  mubeta <- rbind(mu, beta)
  mubeta <- array(mubeta, dim=c(q1, m, ninits))
  
  sigma <- matrix(runif(m*ninits, min=0.01, max=1), nrow=m)*stdR
  
  list(alpha = alpha, mubeta = mubeta, sigma = sigma, gamma = gamma)
  
}  # end function regmixPMLEinit

#' Computes the critical values of the modified EM test.
#' @export
#' @title regmixCrit
#' @name regmixCrit
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma 
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m), 
#' sigma = (sigma_1, ..., sigma_m), gamma = (gamma_1, ..., gamma_m))
#' @param z n by p matrix of regressor associated with gamma
#' @param values 3 by 1 Vector of length 3 (k = 1, 2, 3) at which the p-values are computed
#' @param parallel Determines whether package \code{doParallel} is used for calculation
#' @param cl Cluster used for parallelization; if it is \code{NULL}, the system will automatically
#' @param nrep The number of replications used to compute p-values
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
  if (test.on) # initial values controlled by test.on
    set.seed(test.seed)
  
  y  <- as.vector(y)
  n  <- length(y)
  p  <- 0
  
  alpha   <- parlist$alpha
  mubeta  <- parlist$mubeta
  sigma   <- parlist$sigma
  gamma   <- parlist$gamma
  m       <- length(alpha)
  
  if (!is.null(z)) {
    z <- as.matrix(z)
    p <- ncol(z)
    y   <- y - z %*% gamma
  }
  
  pvals <- NULL
  
  x     <- as.matrix(x)
  x1    <- cbind(1,x)
  q     <- ncol(x)
  
  W  <- t(t(y - x1 %*% mubeta)/sigma)       # normalized data, n by m
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
    S_gamma <- z*rowSums(S_mu)
    S_eta <- cbind(S_eta,S_gamma)
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

#' Computes objective function used in computing LR_2
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

#' Computes Jacobian of the objective function used in computing LR_2
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

#' Computes LR_2 for given Z and I, where q is dim(x)
LR_2.comp <- function(Z, I, q, ninits = 25) {  
  if (test.on) # initial values controlled by test.on
    set.seed(test.seed)
  
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
