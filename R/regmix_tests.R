#' @description Sequentially performs the modified EM test given the data for y and x on 
#' the null hypothesis H_0: m = m_0 where m_0 is in {1, 2, ..., maxm};
#' using this function is equivalent to calling normalmixMEMtestSeq with regressors 
#' specified by x as a parameter.
#' @export
#' @title regmixMEMtestSeq
#' @name regmixMEMtestSeq
#' @param y n by 1 vector of data for y.
#' @param x n by q matrix of data for x.
#' @param z n by p matrix of regressor associated with gamma.
#' @param maxm maximum number of components set as null hypothesis in the mixture.
#' @param ninits number of randomly drawn initial values.
#' @param maxit maximum number of iterations.
#' @param nbtsp number of bootstrap replicates. Default is 199.
#' @param parallel Determines what percentage of available cores are used, 
#' represented by a double in [0,1]. Default is 1.
#' @param cl cluster used for parallelization; if it is \code{NULL}, the system will automatically
#' create a new one for computation accordingly.
#' @param crit.bootstrap.from minimum m in null hypothesis to have critical values 
#' calculated from bootstrap for the test statistics.
#' @param significance.level significance level used for rejecting a null hypothesis.
#' @return A list of with the following items:
#' \item{alpha}{\code{maxm} by \code{maxm} matrix, whose i-th column is a vector of alphas estimated under the null hypothesis m_0 = i.}
#' \item{mu}{\code{maxm} by \code{maxm} matrix, whose i-th column is a vector of mus estimated under the null hypothesis m_0 = i.}
#' \item{sigma}{\code{maxm} by \code{maxm} matrix, whose i-th column is a vector of sigmas estimated under the null hypothesis m_0 = i.}
#' \item{beta}{A list of length \code{maxm}, whose i-th element is a q times i matrix of betas estimated under the null hypothesis m_0 = i.}
#' \item{gam}{\code{maxm} by \code{maxm} matrix, whose i-th column is a vector of gams estimated under the null hypothesis m_0 = i.}
#' \item{emstat}{\code{maxm} by 3 matrix whose i-th row is the value of the modified EM statistic for testing the null hypothesis \eqn{m_0 = i} at k = 1, 2, 3.}
#' \item{pvals}{\code{maxm} by 3 matrix whose i-th row indicates a vector of p-values at k = 1, 2, 3.}
#' \item{aic}{vector of Akaike Information Criterion of the fitted model at m_0 = 1, 2, ..., maxm.}
#' \item{bic}{vector of Bayesian Information Criterion of the fitted model at m_0 = 1, 2, ..., maxm.}
#' \item{loglik}{vector of log-likelihood values of the model at m_0 = 1, 2, ..., maxm.}
#' \item{penloglik}{vector of penalized log-likelihood values of the model at m_0 = 1, 2, ..., maxm.}
#' \item{pmle.result}{list of output from regmixPMLE under the number of components selected by sequantial hypothesis testing.}
#' @examples
#' data(faithful)
#' attach(faithful)
#' \dontrun{regmixMEMtestSeq(y = eruptions, x = waiting, parallel = 1)}
regmixMEMtestSeq <- function (y, x, z = NULL, maxm = 3, ninits = 10, maxit = 2000,
                              nbtsp = 199, parallel = 1, cl = NULL,
                              crit.bootstrap.from = 2,
                              significance.level = 0.05) {
  # Compute the modified EM test statistic for testing H_0 of m components
  # against H_1 of m+1 components for a univariate finite mixture of normals
  if (significance.level >= 1 || significance.level < 0)
    stop ("significance level is not valid.")
  
  y   <- as.vector(y)
  x		<- as.matrix(x)
  n   <- length(y)
  p   <- 0
  q   <- ncol(x)

  if (!is.null(z)) {
    z <- as.matrix(z)
    p <- ncol(z)
    gam <- matrix(0, nrow = p, ncol = maxm)
  }
   else 
    gam <- NULL

  out   <- vector('list', length = maxm)
  aic    <- bic <- double((maxm + 1))
  pvals   <- emstat <- matrix(0, nrow = maxm, ncol = 3)
  loglik  <- penloglik <- double((maxm + 1))

  alpha   <- mu <- sigma <- matrix(0, nrow = maxm, ncol = maxm)
  beta <- list()
  
  # Test H_0:m=1, H_0:m=2, ...
  for (m in 1:maxm){
    pmle.result   <- regmixPMLE(y = y, x = x, m = m, z = z, vcov.method = "none",
                                ninits = ninits, maxit = maxit, binit = binit)
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
    gam0  <- parlist$gam

    alpha[,m] <- c(alpha0, double(maxm - m))
    mu[,m]    <- c(mu0, double(maxm - m))
    sigma[,m] <- c(sigma0, double(maxm - m))
    beta[[m]]  <- c(beta0, double(maxm - m))

    cat(sprintf("%d-component model estimate:\n",m))

    if  (is.null(z)){
      tab = as.matrix(rbind(alpha0, mu0, sigma0, beta0))
      rownames(tab) <- c("alpha", "mu", "sigma", paste("beta", ".", 1:q, sep = ""))
    } else {
      tab = as.matrix(rbind(alpha0, mu0, sigma0, beta0, gam0))
      rownames(tab) <- c("alpha", "mu", "sigma", paste("beta", ".", 1:q, sep = ""), paste("gamma", ".", 1:q, sep = ""))
    }
    colnames(tab) <- c(paste("comp", ".", 1:m, sep = ""))
    print(tab, digits = 4)

    if (!is.null(z)){
      gam[, m] <- gam0
    }
    cat(sprintf("\nAIC, BIC, and log-likelihood of 1 to %.i", m), "component models \n")
    cat(c("AIC    =", sprintf(' %.2f', aic[1:m])), "\n")
    cat(c("BIC    =", sprintf(' %.2f', bic[1:m])), "\n")
    cat(c("loglik =", sprintf('%.2f', loglik[1:m])), "\n\n")

    if (m <= maxm){

      cat(sprintf("Testing the null hypothesis of %d components\n", m))

      an    <- anFormula(parlist = parlist, m = m, n = n)
      par1  <- regmixMaxPhi(y = y, x = x, parlist = parlist, z = z, an = an,
                            ninits = ninits, maxit = maxit)
      emstat.m <- 2*(par1$penloglik-loglik0)
      
      # use the estimate of b as one of the initial values
      binit <- par1$coefficient

      cat(c("modified EM-test statitic ", sprintf('%.3f',emstat.m)),"\n")
      if (m <= crit.bootstrap.from) {
        em.out <- regmixCrit(y = y, x = x, parlist = parlist, z = z,
                             parallel =  parallel, values = emstat.m)
        cat(c("asymptotic p-value       ", sprintf('%.3f',em.out$pvals)),"\n \n")
      } else {
        em.out <- regmixCritBoot(y = y, x = x, parlist=parlist, z = z,
                                 values = emstat.m, ninits = ninits, nbtsp = nbtsp,
                                 parallel = parallel, cl = cl)
        cat(c("bootstrap p-value        ", sprintf('%.3f',em.out$pvals)),"\n \n")
      }
      # noncg.rate[m]   <- par1$noncg.rate
      pvals[m,]     <- em.out$pvals
      emstat[m,]    <- emstat.m
    }
  }
  if (m == maxm)
  {
    pmle.result1   <- regmixPMLE(y = y, x = x, m = (m+1), z = z, vcov.method = "none",
                                    ninits = 2, maxit = maxit, binit = binit)
    loglik[(m+1)] <- pmle.result1$loglik
    penloglik[(m+1)] <- pmle.result1$penloglik
    aic[(m+1)]  <- pmle.result1$aic
    bic[(m+1)]  <- pmle.result1$bic
  }
  
  for (m in 1:maxm) {
    if ( pvals[m,2] >= significance.level ) {
      cat(sprintf("\nThe number of components selected by Sequential Hypothesis Testing (alpha = %.2f) = %.i", 
                  significance.level, m), " \n")     
      cat(sprintf("The number of components selected by AIC = %.i", which.min(aic)), " \n")
      cat(sprintf("The number of components selected by BIC = %.i", which.min(bic)), " \n")
      mubeta <- matrix(0, nrow = (q+1), ncol = m)
      for (j in 1:m) {
        mubeta[1,j] <- mu[j]
        mubeta[-1,j] <- as.vector(beta[[m]][1:q])
      }
      binit <- as.vector(c(alpha[1:m,m], as.vector(mubeta), sigma[1:m,m],  gam[,m]))
      pmle.result   <- regmixPMLE(y = y, x = x, m = m, z = z,
                                     ninits = 2, maxit = maxit, binit = binit)
      cat(sprintf("\nThe summary of the estimated %.i", m), "component model: \n")
      print(summary(pmle.result))
      break
    }
  }
  
  a = list(alpha = alpha, mu = mu, sigma = sigma, beta = beta, gam = gam,
           emstat = emstat, pvals = pvals, aic = aic, bic = bic,
           loglik = loglik, penloglik = penloglik, pmle.result = pmle.result)

  a
}  # end regmixMEMtestSeq

#' @description Performs the modified EM test given the data for y and x on 
#' the null hypothesis H_0: m = m_0. Using this function is equivalent to 
#' calling normalmixMEMtest with regressors specified by x as a parameter.
#' @export
#' @title regmixMEMtest
#' @name regmixMEMtest
#' @param y n by 1 vector of data for y.
#' @param x n by q matrix of data for x.
#' @param m number of components in the mixture defined by a null hypothesis, m_0.
#' @param z n by p matrix of regressor associated with gamma.
#' @param tauset set of initial tau value candidates.
#' @param an tuning parameter used in the penalty function.
#' @param ninits number of randomly drawn initial values.
#' @param crit.method method used to compute the critical values, one of \code{"none"},
#' \code{"asy"}, and \code{"boot"}. The default option is \code{"asy"}. When \code{method = "asy"},
#' the p-values are computed using the asymptotic critical values. When \code{method = "boot"},
#' the p-values are computed by bootstrap.
#' @param nbtsp The number of bootstrap replicates. Default is 199.
#' @param cl cluster used for parallelization; if it is \code{NULL}, the system
#' will automatically generate a new one for computation accordingly.
#' @param parallel Determines what percentage of available cores are used, 
#' represented by a double in [0,1]. Default is 1.
#' @return A list of class \code{normalregMix} with items:
#' \item{emstat}{vector of the value of the modified EM test statistic at k=1, 2, 3.}
#' \item{pvals}{vector of p-values at k = 1, 2, 3.}
#' \item{crit}{When \code{crit.method = "asy"}, \code{crit} is vector of critical values at the 0.1, 0.05, 0.01 level.
#'  When \code{crit.method = "boot"}, \code{crit} is 3 by 3 matrix of critival values at the (0.1, 0.05, 0.01) level, jth row corresponding to k=j.}
#' \item{crit.method}{method used to compute the variance-covariance matrix.}
#' \item{parlist}{parameter estimates as a list containing alpha, mubeta, and sigma (and gam if z is included in the model).}
#' \item{call}{The matched call.}
#' \item{m}{number of components in the mixture defined by the null hypothesis, \eqn{m_0}.}
#' @examples
#' data(faithful)
#' attach(faithful)
#' \dontrun{regmixMEMtest(y = eruptions, x = waiting, m = 1, crit.method = "asy")}
#' \dontrun{regmixMEMtest(y = eruptions, x = waiting, m = 2, crit.method = "asy")}
regmixMEMtest <- function (y, x, m = 2, z = NULL, tauset = c(0.1,0.3,0.5),
                           an = NULL, ninits = 100,
                           crit.method = c("none", "asy", "boot"), nbtsp = 199,
                           cl = NULL,
                           parallel = 1) {
  # Compute the modified EM test statistic for testing H_0 of m components
  # against H_1 of m+1 components for a univariate finite mixture of normals
  y <- as.vector(y)
  x <- as.matrix(x)
  q <- ncol(x)
  n <- length(y)

  if (!is.null(z))
    z <- as.matrix(z)

  regmix.pmle.result    <- regmixPMLE(y=y, x=x, m=m, z=z, vcov.method="none", ninits=ninits)
  loglik0 <- regmix.pmle.result$loglik

  if (is.null(an))
    an <- anFormula(parlist = regmix.pmle.result$parlist, m = m, n = n, q = q)

  par1    <- regmixMaxPhi(y=y, x=x, parlist=regmix.pmle.result$parlist, z=z,
                          an=an, tauset = tauset, ninits=ninits)
  
  
  emstat <- 2*(par1$penloglik-loglik0)

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
