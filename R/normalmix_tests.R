#' @description  Sequentially performs the modified EM test given the data for y and x 
#' on the null hypothesis \eqn{H_0: m = m_0} where \eqn{m_0} is in (1, 2, ..., maxm)
#' @export
#' @title normalmixMEMtestSeq
#' @name normalmixMEMtestSeq
#' @param y n by 1 vector of data for y.
#' @param x n by q matrix of data for x (if exists).
#' @param z n by p matrix of regressor associated with gamma.
#' @param maxm The maximum number of components in the mixture set as the null hypothesis.
#' @param ninits The number of randomly drawn initial values.
#' @param maxit The maximum number of iterations.
#' @param nbtsp The number of bootstrap observations. Defalt value is 199.
#' @param parallel Determines what percentage of available cores are used, 
#' represented by a double in [0,1]. Default value is 1.
#' @param cl Cluster used for parallelization; if it is \code{NULL}, the system will automatically
#' create a new one for computation accordingly.
#' @param crit.bootstrap.from The minimum m in null hypothesis to have critical values 
#' calculated from bootstrap for the test statistics.
#' @param significance.level Significance level used for rejecting a null hypothesis.
#' @param LRT.penalized Determines whether penalized likelihood is used in calculation of LRT
#' statistic for likelihood in an alternative hypothesis.
#' @return A list of with the following items:
#' \item{alpha}{A \code{maxm} by \code{maxm} matrix, whose i-th column is a vector of alphas estimated under the null hypothesis \eqn{m_0 = i}.}
#' \item{mu}{A \code{maxm} by maxm \code{maxm}, whose i-th column is a vector of mus estimated under the null hypothesis \eqn{m_0 = i}.}
#' \item{sigma}{A \code{maxm} by \code{maxm} matrix, whose i-th column is a vector of sigmas estimated under the null hypothesis \eqn{m_0 = i}.}
#' \item{beta}{A list of length \code{maxm}, whose i-th element is a q times i matrix of betas estimated under the null hypothesis \eqn{m_0 = i}.}
#' \item{gam}{A p by \code{maxm} matrix, whose i-th column is gamma estimated under the null hypothesis \eqn{m_0 = i}.}
#' \item{emstat}{A \code{maxm} by 3 matrix whose i-th row is the value of the modified EM statistic for testing the null hypothesis \eqn{m_0 = i} at k = 1, 2, 3.}
#' \item{pvals}{A \code{maxm} by 3 matrix whose i-th row indicates a vector of p-values at k = 1, 2, 3.}
#' \item{aic}{The vector of the values of Akaike Information Criterion of the fitted model at \eqn{m_0 = 1, 2, ..., maxm}.}
#' \item{bic}{The vector of the values of Bayesian Information Criterion of the fitted model at \eqn{m_0 = 1, 2, ..., maxm}.}
#' \item{loglik}{The vector of log-likelihood values of the fitted model at \eqn{m_0 = 1, 2, ..., maxm}.}
#' \item{penloglik}{The vector of penalized log-likelihood values of the fitted model at \eqn{m_0 = 1, 2, ..., maxm}.}
#' \item{pmle.result}{A list of output from normalmixPMLE under the number of components selected by sequantial hypothesis testing.}
#' @examples
#' data(faithful)
#' attach(faithful)
#' normalmixMEMtestSeq(y = eruptions, parallel = 0)
normalmixMEMtestSeq <- function (y, x = NULL, z = NULL,  maxm = 3, ninits = 10, maxit = 2000,
                                 nbtsp = 199, parallel = 1, cl = NULL,
                                 crit.bootstrap.from = 3, significance.level = 0.05,
                                 LRT.penalized = FALSE) {
  # Compute the modified EM test statistic for testing H_0 of m components
  # against H_1 of m+1 components for a univariate finite mixture of normals

  if (significance.level >= 1 || significance.level < 0)
    stop ("significance level is not valid.")
  if (!is.null(x))
    return (regmixMEMtestSeq(y = y, x = x, z = z, maxm = maxm, ninits = ninits, maxit = maxit,
                             nbtsp = nbtsp, parallel = parallel, cl = cl,
                             crit.bootstrap.from = crit.bootstrap.from))
  
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

      an    <- anFormula(parlist = parlist, m = m, n = n, LRT.penalized = LRT.penalized)
      par1  <- normalmixMaxPhi(y = y, parlist = parlist, z = z, an = an,
                               ninits = ninits, maxit = maxit)
      emstat.m  <- 2*(par1$loglik - loglik0)
      if (LRT.penalized) # use the penalized log-likelihood.
        emstat.m  <- 2*(par1$penloglik - loglik0)

      # use the estimate of b as one of the initial values
      binit <- par1$coefficient

      cat(c("modified EM-test statitic ", sprintf('%.3f',emstat.m)),"\n")
      if (m <= crit.bootstrap.from) {
        em.out <- normalmixCrit(y=y, parlist=parlist, z=z, values = emstat.m)
        cat(c("asymptotic p-value       ", sprintf('%.3f',em.out$pvals)),"\n \n")
      } else {
        em.out <- normalmixCritBoot(y=y, parlist=parlist, z=z, values = emstat.m,
                                    ninits = ninits, nbtsp = nbtsp, parallel = parallel, cl = cl,
                                    LRT.penalized = LRT.penalized)
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

#' Performs the modified EM test given the data for y and x on the null hypothesis \eqn{H_0: m = m_0}.
#' @export
#' @title normalmixMEMtest
#' @name normalmixMEMtest
#' @param y n by 1 vector of data.
#' @param x n by q matrix of data for x (if exists).
#' @param m The number of components in the mixture defined by the null hypothesis, \eqn{m_0}.
#' @param z n by p matrix of regressor associated with gamma.
#' @param tauset A set of initial tau value candidates.
#' @param an a term used for penalty function.
#' @param ninits The number of randomly drawn initial values.
#' @param crit.method Method used to compute the variance-covariance matrix, one of \code{"none"},
#' \code{"asy"}, and \code{"boot"}. The default option is \code{"asy"}. When \code{method = "asy"},
#' the p-values are computed based on an asymptotic method. When \code{method = "boot"},
#' the p-values are generated by bootstrapping.
#' @param nbtsp The number of bootstrap observations. Default value is 199.
#' @param cl Cluster used for parallelization; if it is \code{NULL}, the system will automatically
#' create a new one for computation accordingly.
#' @param parallel Determines what percentage of available cores are used, 
#' represented by a double in [0,1]. Default value is 1.
#' @param LRT.penalized Determines whether penalized likelihood is used in calculation of LRT
#' statistic for likelihood in an alternative hypothesis.
#' @return A list of class \code{normalMix} with items:
#' \item{coefficients}{A vector of parameter estimates. Ordered as \eqn{\alpha_1,\ldots,\alpha_m,\mu_1,\ldots,\mu_m,\sigma_1,\ldots,\sigma_m,gam}, where gam is the coefficient of z.}
#' \item{parlist}{The parameter estimates as a list containing alpha, mu, and sigma (and gam if z is included in the model).}
#' \item{loglik}{The maximized value of the log-likelihood.}
#' \item{penloglik}{The maximized value of the penalized log-likelihood.}
#' \item{aic}{Akaike Information Criterion of the fitted model.}
#' \item{bic}{Bayesian Information Criterion of the fitted model.}
#' \item{postprobs}{n by m matrix of posterior probabilities for observations.}
#' \item{components}{n by 1 vector of integers that indicates the indices of components
#' each observation belongs to based on computed posterior probabilities.}
#' \item{call}{The matched call.}
#' \item{m}{The number of components in the mixture defined by the null hypothesis, \eqn{m_0}.}
#' @examples
#' data(faithful)
#' attach(faithful)
#' normalmixMEMtest(y = eruptions, m = 1, crit.method = "asy", parallel = 0)
#' normalmixMEMtest(y = eruptions, m = 2, crit.method = "asy", parallel = 0)
normalmixMEMtest <- function (y, x = NULL, m = 2, z = NULL, an = NULL, tauset = c(0.1,0.3,0.5),
                              ninits = 10,
                              crit.method = c("asy", "boot", "none"), nbtsp = 199,
                              cl = NULL,
                              parallel = 1,
                              LRT.penalized = FALSE) {
  # Compute the modified EM test statistic for testing H_0 of m components
  # against H_1 of m+1 components for a univariate finite mixture of normals
  if (!is.null(x))
    return (regmixMEMtest(y = y, x = x, m = m, z = z, tauset = tauset, an = an,
                          ninits = ninits, crit.method = crit.method, nbtsp = nbtsp,
                          cl = cl, parallel = parallel))
  y     <- as.vector(y)
  n     <- length(y)
  if (!is.null(z)) {z <- as.matrix(z)}
  crit.method <- match.arg(crit.method)

  pmle.result    <- normalmixPMLE(y=y, m=m, z=z, vcov.method = "none", ninits=ninits)
  loglik0 <- pmle.result$loglik

  if (is.null(an)){ an <- anFormula(parlist = pmle.result$parlist, m = m, n = n, 
                                    LRT.penalized = LRT.penalized) }

  par1  <- normalmixMaxPhi(y=y, parlist=pmle.result$parlist, z=z,
                           an=an, tauset = tauset, ninits=ninits)

  emstat  <- 2*(par1$loglik - loglik0)
  if (LRT.penalized) # use the penalized log-likelihood.
    emstat  <- 2*(par1$penloglik - loglik0)

  if (crit.method == "asy"){
    result  <- normalmixCrit(y=y, parlist=pmle.result$parlist, z=z, values=emstat)
  } else if (crit.method == "boot") {
    result  <- normalmixCritBoot(y=y, parlist= pmle.result$parlist, z=z, values=emstat,
                                 ninits=ninits, nbtsp=nbtsp, parallel=parallel, cl=cl,
                                 LRT.penalized = LRT.penalized)
  } else {
    result <- list()
    result$crit <- result$pvals <- rep(NA,3)
  }

  a <- list(emstat = emstat, pvals = result$pvals, crit = result$crit, crit.method = crit.method,
            parlist = pmle.result$parlist, call = match.call(), m = m, label = "MEMtest")

  class(a) <- "normalregMix"

  a

}  # end normalmixMEMtest
