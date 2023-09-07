#' @description  Sequentially performs the modified EM test given the data for y 
#' on the null hypothesis \eqn{H_0: m = m_0} where \eqn{m_0} is in (1, 2, ..., maxm)
#' @export
#' @title normalmixMEMtestSeq
#' @name normalmixMEMtestSeq
#' @param y n by 1 vector of data for y.
#' @param z n by p matrix of regressor associated with gamma.
#' @param maxm maximum number of components in the mixture set as the null hypothesis.
#' @param ninits number of randomly drawn initial values.
#' @param maxit maximum number of iterations.
#' @param nbtsp number of bootstrap replicates. Defalt is 199.
#' @param parallel Determines what percentage of available cores are used, represented by a double in [0,1]. Default is 1.
#' @param cl cluster used for parallelization; if it is \code{NULL}, the system will automatically
#' create a new one for computation accordingly.
#' @param crit.bootstrap.from minimum m in null hypothesis to have critical values 
#' calculated from bootstrap for the test statistics.
#' @param significance.level significance level used for rejecting a null hypothesis.
#' @return A list with the following items:
#' \item{alpha}{\code{maxm} by \code{maxm} matrix, whose i-th column is a vector of alphas estimated under the null hypothesis \eqn{m_0 = i}.}
#' \item{mu}{\code{maxm} by maxm \code{maxm}, whose i-th column is a vector of mus estimated under the null hypothesis \eqn{m_0 = i}.}
#' \item{sigma}{\code{maxm} by \code{maxm} matrix, whose i-th column is a vector of sigmas estimated under the null hypothesis \eqn{m_0 = i}.}
#' \item{beta}{A list of length \code{maxm}, whose i-th element is a q times i matrix of betas estimated under the null hypothesis \eqn{m_0 = i}.}
#' \item{gam}{p by \code{maxm} matrix, whose i-th column is gamma estimated under the null hypothesis \eqn{m_0 = i}.}
#' \item{emstat}{\code{maxm} by 3 matrix whose i-th row is the value of the modified EM statistic for testing the null hypothesis \eqn{m_0 = i} at k = 1, 2, 3.}
#' \item{pvals}{\code{maxm} by 3 matrix whose i-th row indicates a vector of p-values at k = 1, 2, 3.}
#' \item{aic}{vector of the values of Akaike Information Criterion of the fitted model at \eqn{m_0 = 1, 2, ..., maxm}.}
#' \item{bic}{vector of the values of Bayesian Information Criterion of the fitted model at \eqn{m_0 = 1, 2, ..., maxm}.}
#' \item{loglik}{vector of log-likelihood values of the fitted model at \eqn{m_0 = 1, 2, ..., maxm}.}
#' \item{penloglik}{vector of penalized log-likelihood values of the fitted model at \eqn{m_0 = 1, 2, ..., maxm}.}
#' \item{pmle.result}{list of output from normalmixPMLE under the number of components selected by sequantial hypothesis testing.}
#' @examples
#' data(faithful)
#' attach(faithful)
#' normalmixMEMtestSeq(y = eruptions, parallel = 0)
normalmixMEMtestSeq <- function (y, z = NULL, maxm = 3, ninits = 10, maxit = 2000,
                                 nbtsp = 199, parallel = 1, cl = NULL,
                                 crit.bootstrap.from = 3, significance.level = 0.05) {
  # Compute the modified EM test statistic for testing H_0 of m components
  # against H_1 of m+1 components for a univariate finite mixture of normals

  if (significance.level >= 1 || significance.level < 0)
    stop ("significance level is not valid.")
  
  y   <- as.vector(y)
  n   <- length(y)
  p   <- 0
  if (!is.null(z)) {
    z <- as.matrix(z)
    p <- ncol(z)
    gam <- matrix(0, nrow = p, ncol = maxm)
  }
  else 
    gam <- NULL
	
  out   <- vector('list', length = maxm)
  aic    <- bic <- double((maxm+1))
  pvals   <- emstat <- matrix(0, nrow = maxm, ncol = 3)
  loglik  <- penloglik <- double((maxm + 1))

  alpha   <- mu <- sigma <- matrix(0, nrow = maxm, ncol = maxm)

  # Test H_0:m=1, H_0:m=2, ...
  binit <- NULL
  for (m in 1:maxm){
    pmle.result   <- normalmixPMLE(y = y, m = m, z = z, vcov.method = "none",
                                   ninits = ninits, maxit = maxit, binit = binit)
    loglik[m] <- loglik0 <- pmle.result$loglik
    penloglik[m] <- penloglik0 <- pmle.result$penloglik
    aic[m]  <- pmle.result$aic
    bic[m]  <- pmle.result$bic

    parlist <- pmle.result$parlist
    alpha0  <- parlist$alpha
    mu0     <- parlist$mu
    sigma0  <- parlist$sigma
    gam0  <- parlist$gam

    alpha[,m] <- c(alpha0, double(maxm - m))
    mu[,m]     <- c(mu0, double(maxm - m))
    sigma[,m] <- c(sigma0, double(maxm - m))

    cat(sprintf("%d-component model estimate:\n",m))
    tab = as.matrix(rbind(alpha0, mu0, sigma0))
    rownames(tab) <- c("alpha", "mu", "sigma")
    colnames(tab) <- c(paste("comp", ".", 1:m, sep = ""))
    print(tab, digits = 4)

    if (!is.null(z)){
      gam[, m] <- gam0
      cat("gam =", gam0,"\n")
    }
    cat(sprintf("\nAIC, BIC, and log-likelihood of 1 to %.i", m), "component models \n")
    cat(c("AIC    =", sprintf(' %.2f', aic[1:m])), "\n")
    cat(c("BIC    =", sprintf(' %.2f', bic[1:m])), "\n")
    cat(c("loglik =", sprintf('%.2f', loglik[1:m])), "\n\n")

    if (m <= maxm){

      cat(sprintf("Testing the null hypothesis of %d components\n", m))

      an    <- anFormula(parlist = parlist, m = m, n = n)
      cn <- min(m,4)
      par1  <- normalmixMaxPhi(y = y, parlist = parlist, z = z, an = an, cn = cn,
                               ninits = ninits, maxit = maxit)
      emstat.m  <- 2*(par1$penloglik - loglik0)

      # use the estimate of b as one of the initial values
      binit <- par1$coefficient

      cat(c("modified EM-test statitic ", sprintf('%.3f',emstat.m)),"\n")
      if (m <= crit.bootstrap.from) {
        em.out <- normalmixCrit(y=y, parlist=parlist, z=z, values = emstat.m)
        cat(c("asymptotic p-value       ", sprintf('%.3f',em.out$pvals)),"\n \n")
      } else {
        em.out <- normalmixCritBoot(y=y, parlist=parlist, z=z, values = emstat.m,
                                    ninits = ninits, nbtsp = nbtsp, parallel = parallel, cl = cl)
        cat(c("bootstrap p-value        ", sprintf('%.3f',em.out$pvals)),"\n \n")
      }
      pvals[m,]     <- em.out$pvals
      emstat[m,]    <- emstat.m
    }
    if (m == maxm)
    {
      pmle.result1   <- normalmixPMLE(y = y, m = (m+1), z = z, vcov.method = "none",
                                     ninits = 2, maxit = maxit, binit = binit)
      loglik[(m+1)] <- pmle.result1$loglik
      penloglik[(m+1)] <- pmle.result1$penloglik
      aic[(m+1)]  <- pmle.result1$aic
      bic[(m+1)]  <- pmle.result1$bic
    }

  }

  for (m in 1:maxm)
    if ( pvals[m,2] >= significance.level ) {
      cat(sprintf("\nThe number of components selected by Sequential Hypothesis Testing (alpha = %.2f) = %.i", 
                  significance.level, m), " \n")
      cat(sprintf("The number of components selected by AIC = %.i", which.min(aic)), " \n")
      cat(sprintf("The number of components selected by BIC = %.i", which.min(bic)), " \n")
      binit <- as.vector(c(alpha[1:m,m], mu[1:m,m], sigma[1:m,m],  gam[,m]))
      pmle.result   <- normalmixPMLE(y = y, m = m, z = z,
                                     ninits = 2, maxit = maxit, binit = binit)
      cat(sprintf("\nThe summary of the estimated %.i", m), "component model: \n")
      print(summary(pmle.result))
      break
    }

  a = list(alpha = alpha, mu = mu, sigma = sigma, gam = gam, emstat = emstat, pvals = pvals, aic = aic, bic = bic, loglik = loglik, penloglik = penloglik, pmle.result = pmle.result)

  a
}  # end normalmixMEMtestSeq


#' Performs the modified EM test given the data for y on the null hypothesis \eqn{H_0: m = m_0}.
#' @export
#' @title normalmixMEMtest
#' @name normalmixMEMtest
#' @param y n by 1 vector of data.
#' @param m number of components in the mixture defined by the null hypothesis, \eqn{m_0}.
#' @param z n by p matrix of regressor associated with gamma.
#' @param an tuning parameter used in the penalty function on sigma.
#' @param ninits number of randomly drawn initial values.
#' @param crit.method method used to compute the critical values, one of \code{"none"},
#' \code{"asy"}, and \code{"boot"}. The default option is \code{"asy"}. When \code{method = "asy"}
#' and \eqn{m_0 \le 3}, the p-values are computed using the asymptotic critical values.
#' When \code{method = "boot"} or \eqn{m_0 \ge 4}, the p-values are computed by bootstrap.
#' @param nbtsp number of bootstrap replicates. Default value is 199.
#' @param cl cluster used for parallelization; if \code{NULL}, the system will automatically
#' create a new one for computation accordingly.
#' @param parallel Determines the percentage of available cores used, represented by a double in [0,1]. Default is 1.
#' @return A list of class \code{normalMix} with items:
#' \item{emstat}{vector of the value of the modified EM test statistic at k=1, 2, 3.}
#' \item{pvals}{vector of p-values at k = 1, 2, 3.}
#' \item{crit}{When \code{crit.method = "asy"}, \code{crit} is vector of critical values at the 0.1, 0.05, 0.01 level.
#'  When \code{crit.method = "boot"}, \code{crit} is 3 by 3 matrix of critival values at the (0.1, 0.05, 0.01) level, jth row corresponding to k=j.}
#' \item{crit.method}{method used to compute the variance-covariance matrix.}
#' \item{parlist}{parameter estimates as a list containing alpha, mu, and sigma (and gam if z is included in the model).}
#' \item{call}{The matched call.}
#' \item{m}{number of components in the mixture defined by the null hypothesis, \eqn{m_0}.}
#' @examples
#' data(faithful)
#' attach(faithful)
#' normalmixMEMtest(y = waiting, m = 1, crit.method = "asy", parallel = 0)
#' normalmixMEMtest(y = waiting, m = 2, crit.method = "asy", parallel = 0)
normalmixMEMtest <- function (y, m = 2, z = NULL, an = NULL, ninits = 10,
                              crit.method = c("asy", "boot", "none"), nbtsp = 199,
                              cl = NULL, parallel = 1) {
  # Compute the modified EM test statistic for testing H_0 of m components
  # against H_1 of m+1 components for a univariate finite mixture of normals
  y     <- as.vector(y)
  n     <- length(y)
  if (!is.null(z)) {
    z     <- as.matrix(z)
    p     <- ncol(z)
  } else {
    p <- 0
  }
  
  crit.method <- match.arg(crit.method)
  
  pmle.result    <- normalmixPMLE(y=y, m=m, z=z, vcov.method = "none", ninits=ninits)
  loglik0 <- pmle.result$loglik
  
  if (is.null(an)){ an <- anFormula(parlist = pmle.result$parlist, m = m, n = n) }
  
  cn <- min(m,4)
  tauset <- c(0.1, 0.3, 0.5)
  htau <- t(as.matrix(expand.grid(c(1:m), tauset)))
  
  ninits.maxphi <- floor(ninits/2)
  
  results <- apply(htau, 2, normalmixMaxPhi, y, parlist=pmle.result$parlist,
                   z, an, cn, ninits=ninits.maxphi)
  loglik.all <- t(sapply(results, "[[", "loglik")) # each row = each htau, each col = each k
  penloglik.all <- t(sapply(results, "[[", "penloglik")) # each row = each htau, each col = each k
  coefficient.all <- t(sapply(results, "[[", "coefficient")) # each row = each htau, each col = each k
  
  loglik <- apply(loglik.all, 2, max)  # 3 by 1 vector corresponding to k=1,2,3
  penloglik <- apply(penloglik.all, 2, max)  # 3 by 1 vector corresponding to k=1,2,3
  index <- which.max(loglik.all[ ,3]) # a pair (h,m) that gives the highest likelihood at k=3
  coefficient <- as.vector(coefficient.all[index,])
  
  emstat  <- 2*(penloglik - loglik0)
  
  if (crit.method == "none") {
    result <- list()
    result$crit <- result$pvals <- rep(NA,3)
  } else if ((crit.method == "asy") && (m <=3)) {
    result  <- normalmixCrit(y=y, parlist=pmle.result$parlist, z=z, values=emstat)
  } else {
    crit.method <- "boot"
    result  <- normalmixCritBoot(y=y, parlist= pmle.result$parlist, z=z, values=emstat,
                                 ninits=ninits, nbtsp=nbtsp, parallel=parallel, cl=cl,
                                 an = an)
  }
  
  a <- list(emstat = emstat, pvals = result$pvals, crit = result$crit, crit.method = crit.method,
            parlist = pmle.result$parlist, call = match.call(), m = m, label = "MEMtest")
  
  class(a) <- "normalregMix"
  
  a
  
}  # end normalmixMEMtest

#' @description Given a pair of h and tau and data, compute ordinary &
#' penalized log-likelihood ratio resulting from MEM algorithm at k=1,2,3, 
#' tailored for parallelization.
#' @title normalmixMaxPhi
#' @name normalmixMaxPhi
#' @param htaupair set of h and tau.
#' @param y n by 1 vector of data.
#' @param parlist parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam).
#' @param z n by p matrix of regressor associated with gamma.
#' @param an tuning parameter used in the penalty function on sigma.
#' @param cn tuning parameter used in the penalty function on tau.
#' @param ninits number of randomly drawn initial values.
#' @param epsilon.short convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param epsilon convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}.
#' @param maxit.short maximum number of iterations in short EM.
#' @param maxit maximum number of iterations.
#' @param verb Determines whether to print a message if an error occurs.
#' @return A list of phi, log-likelihood, and penalized log-likelihood resulting from MEM algorithm.
normalmixMaxPhi <- function (htaupair, y, parlist, z = NULL, an, cn,
                             ninits, epsilon.short = 1e-02, epsilon = 1e-08,
                             maxit.short = 500, maxit = 2000, verb = FALSE)
{
  alpha0 <- parlist$alpha
  
  m      <- length(alpha0)
  m1     <- m+1
  k      <- 1
  n      <- length(y)
  h      <- as.numeric(htaupair[1])
  tau    <- as.numeric(htaupair[2])
  
  mu0    <- parlist$mu
  mu0h   <- c(-1e+10,mu0,1e+10)        # m+2 by 1
  sigma0 <- parlist$sigma
  sigma0h<- c(sigma0[1:h],sigma0[h:m]) # m+1 by 1
  gam0 <- parlist$gam
  if (is.null(z)) {
    ztilde <- matrix(0) # dummy
    gam <- NULL
    p <- 0
  } else {
    ztilde <- as.matrix(z)
    p      <- ncol(z)
  }
  
  ninits.short <- ninits*10*(1+p)*m
  
  # Generate initial values
  tmp <- normalmixPhiInit(y = y, parlist = parlist, z = z, h=h, tau = tau, ninits = ninits.short)
  
  # short EM
  b0 <- as.matrix(rbind(tmp$alpha, tmp$mu, tmp$sigma, tmp$gam))
  out.short <- cppNormalmixPMLE(b0, y, ztilde, mu0h, sigma0h, m1, p, an, cn, maxit.short, ninits.short,
                                epsilon.short, tau, h, k)
  components <- order(out.short$penloglikset, decreasing = TRUE)[1:ninits]
  if (verb && any(out.short$notcg)) {
    cat(sprintf("non-convergence rate at short-EM = %.3f\n",mean(out.short$notcg)))
  }
  # long EM
  b1 <- as.matrix(b0[ ,components])
  out <- cppNormalmixPMLE(b1, y, ztilde, mu0h, sigma0h, m1, p, an, cn, maxit, ninits, epsilon, tau, h, k)
  
  index     <- which.max(out$penloglikset)
  alpha <- b1[1:m1,index]
  mu <- b1[(1+m1):(2*m1),index]
  sigma <- b1[(1+2*m1):(3*m1),index]
  if (!is.null(z)) {
    gam     <- b1[(3*m1+1):(3*m1+p),index]
  }
  mu.order  <- order(mu)
  alpha     <- alpha[mu.order]
  mu        <- mu[mu.order]
  sigma     <- sigma[mu.order]
  sigma0h <- sigma0h[mu.order]
  b <- as.matrix( c(alpha, mu, sigma, gam) )
  
  # initilization for taking two EM steps
  loglik <-  vector("double", 3)
  penloglik <-  vector("double", 3)
  coefficient <- vector("double", length(b))
  
  penloglik[1] <- out$penloglikset[[index]]
  loglik[1]    <- out$loglikset[[index]]
  for (k in 2:3) {
    ninits <- 1
    maxit <- 1  # SHOULD BE 1
    # Two EM steps
    out <- cppNormalmixPMLE(b, y, ztilde, mu0h, sigma0h, m1, p, an, cn, maxit, ninits, epsilon, tau, h, k)
    alpha <- b[1:m1,1] # b has been updated
    mu <- b[(1+m1):(2*m1),1]
    sigma <- b[(1+2*m1):(3*m1),1]
    if (!is.null(z)) {
      gam     <- b[(3*m1+1):(3*m1+p),1]
    }
    loglik[k]     <- out$loglikset[[1]]
    penloglik[k]  <- out$penloglikset[[1]]
    
    # Check singularity: if singular, break from the loop
    if ( any(sigma < 1e-06) || any(alpha < 1e-06) || is.na(sum(alpha)) ) {
      loglik[k]    <- -Inf
      penloglik[k] <- -Inf
      break
    }
    
    mu.order  <- order(mu)
    alpha     <- alpha[mu.order]
    mu        <- mu[mu.order]
    sigma     <- sigma[mu.order]
    sigma0h <- sigma0h[mu.order]
  }
  coefficient <- as.matrix( c(alpha, mu, sigma, gam) ) # at k=3
  
  return (list(coefficient = coefficient, loglik = loglik, penloglik = penloglik))
}  # end function normalmixMaxPhi
