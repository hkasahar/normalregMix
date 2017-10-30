#' @description Turns on/off the test mode.
#' @export
#' @title normalregMixTestMode
#' @name normalregMixTestMode
#' @description \code{normalregMix} uses random numbers
#' (i) to draw the intial values of the penalized MLE and the
#' modified EM test, and (ii) to draw bootstrapped samples in bootstrapped tests.
#' When the test mode is turned on, these random numbers are drawn 
#' with the random seed provided. This method would be suitable
#' for users who would like to replicate experiments. By default, the test mode is turned off.
#' @param on option to turn on the test mode.
#' @param seed random seed to be used for initialization.
#' @param hide.message Determines whether to print the current seed and status.
#' @examples
#' normalregMixTestMode(TRUE)
#' normalregMixTestMode(FALSE)
normalregMixTestMode <- function(on = FALSE, seed = 8888577, hide.message = FALSE)
{
  normalregMixtest.env$normalregMix.test.on <- on
  normalregMixtest.env$normalregMix.test.seed <- seed
  
  if (!hide.message)
    print(paste("The test mode is currently",
                switch(as.character(normalregMixtest.env$normalregMix.test.on), "TRUE" = "ON", "FALSE" = "OFF"),
                "with normalregMixtest.env$normalregMix.test.seed",
                as.character(normalregMixtest.env$normalregMix.test.seed)))
}


#' @description Generates a vector that indicates which component each observation belongs to,
#' based on its posterior probability.
#' @export
#' @title getComponentcomponents
#' @name getComponentcomponents
#' @param postprobs n by m matrix of posterior probabilities for
#' m-component model on n observations.
#' @return n by 1 vector of components that indicate which component each
#'  observation belongs to based on its posterior probability.
getComponentcomponents <- function(postprobs)
{
  postprobs.mat <- as.matrix(postprobs)
  apply(postprobs.mat, 1, function(i) (which(i==max(i))))
}

#' @description Returns the summary of a normalregMix instance.
#' @export
#' @title summary.normalregMix
#' @name summary.normalregMix
#' @param object normalregMix instance.
#' @param reorder Determines whether components are reordered in summary.
#' @param digits digits used for reporting.
#' @param ... other arguments that do not affect the function.
summary.normalregMix <- function(object, reorder = FALSE, digits = 3, ...) {

if (object$label == "PMLE") {
  coef <- object$coefficients
  vcov <- object$vcov
  m    <- object$m
  # When m=1, the first element of coef is alpha (=1), so we drop it
  if (object$m == 1) { coef <- coef[-1] }

  if (reorder && m != 1){
    len1 <- length(object$coefficients)
    k <- nrow(object$parlist$mubeta)
    if (is.null(k)) { k <- 1 }
    sel.vec <- NULL
    for (j in 1:m){
      sel <- c(j, (m+(j-1)*k+1):(m+j*k), j+m*(k+1))
      sel.vec <- c(sel.vec,sel)
    }
    if (!is.null(object$parlist$gam)) {
      p <- length(object$parlist$gam)
      sel <- c((len1-p+1):len1)
      sel.vec <- c(sel.vec,sel)
    }
    reorder.mat <- matrix(0, nrow=len1, ncol=len1)
    sel.vec.2 <- cbind(1:len1, sel.vec)
    reorder.mat[sel.vec.2] <- 1
    coef <- coef[sel.vec]
    vcov <- reorder.mat %*% vcov %*% t(reorder.mat)
  }

  se    <- sqrt(diag(vcov))
  tval  <- coef / se
  TAB   <- cbind(Estimate = coef, StdErr =se, t.value = tval, p.value = 2*pnorm(-abs(tval)))
  res   <- list(coefficients = TAB, parlist = object$parlist, vcov = vcov,
              loglik = object$loglik, penloglik = object$penloglik,
              aic = object$aic, bic = object$bic, call = object$call,
              m = object$m, reorder = reorder,
              digits = digits)
  class(res) <- "summary.normalregMix"
  res

} else {
  stop("The object is not an output of normalmixPMLE/regmixPMLE.")
}

}  # end function summary.normalregMix

#' @description Prints the summary of a normalregMix instance.
#' @export
#' @title print.summary.normalregMix
#' @name print.summary.normalregMix
#' @param x normalregMix instance.
#' @param ... Other arguments that do not affect the function.
print.summary.normalregMix <- function(x, ...) {

cat("\nCall:\n")
print(x$call)
cat("\n")

m   <- x$m
tab <- x$coefficients
reorder <- x$reorder
k <- nrow(x$parlist$mubeta)
if (is.null(k)) { k <- 1 }

if (reorder) {
#   cat(sprintf("Number of components: %d\n",m))
  for (j in 1:m){
    cat(sprintf("Component %i\n",j))
    coef.j = tab[c(((k+2)*(j-1)+1):((k+2)*j)), ]
    # rownames(coef.j) <- c("alpha","mu","sigma")
    printCoefmat(coef.j, P.values = TRUE, has.Pvalue = TRUE, signif.legend = FALSE, digits=x$digits)
  }
  if (!is.null(x$parlist$gam)) {
    p <- length(x$parlist$gam)
    gam <- tab[(nrow(tab)-p+1):nrow(tab), , drop = FALSE]
    cat("gam\n")
    printCoefmat(gam, P.values = TRUE, has.Pvalue = TRUE, signif.legend = FALSE, digits=x$digits)
  }
  cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
} else {
  printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE, digits = x$digits)
}

cat(sprintf("\nloglik at estimate:  %.3f\n", x$loglik))
cat(sprintf("penloglik at estimate: %.3f\n", x$penloglik))
cat(sprintf("AIC: %.3f\n", x$aic))
cat(sprintf("BIC: %.3f\n", x$bic))

}

#' @description Prints the description of a normalregMix instance.
#' @export
#' @title print.normalregMix
#' @name print.normalregMix
#' @param x normalregMix instance.
#' @param ... other arguments that do not affect the function.
print.normalregMix <- function(x, ...) {

cat("\nCall:\n")
print(x$call)

if (x$label == "MEMtest") {
  cat(sprintf("\nTesting the null hypothesis of %d components\n", x$m))
  cat("                            k = 1  k = 2  k = 3 \n")
  cat(c("modified EM-test statistic ",sprintf('%.3f ', x$emstat)),"\n")
  if (x$crit.method == "asy") {
  cat(c("asymptotic p-value         ",sprintf('%.3f ', x$pvals)),"\n")
  } else if (x$crit.method == "boot") {
    cat(c("bootstrap p-value          ",sprintf('%.3f ', x$pvals)),"\n")
  }
} else if (x$label == "PMLE") {
  cat("\nCoefficients:\n")
  print(x$coefficients, digits=4)
  cat(sprintf("\nloglik at estimate: %.3f\n", x$loglik))
  cat(sprintf("penloglik at estimate: %.3f\n", x$penloglik))
  cat(sprintf("AIC: %.3f\n", x$aic))
  cat(sprintf("BIC: %.3f\n", x$bic))
} else {
  stop("The object is not a valid normalMix object.")
}
}


#' @description Computes the tuning parameter \eqn{a_n} based on
#' empirical formulas obtained by a similar method to Kasahara and Shimotsu (2015).
#' @export
#' @title anFormula
#' @name anFormula
#' @param parlist parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gamma_1, ..., gamma_m)).
#' @param m number of components in the mixture.
#' @param n number of observations.
#' @param q dimension of x (default is 0).
#' @return tuning parameter \eqn{a_n}.
#' @references Kasahara, H., and Shimotsu, K. (2015)
#' Testing the Number of Components in Normal Mixture Regression Models,
#' \emph{Journal of the American Statistical Association},
#' \bold{110}, 1632--1645.
anFormula <- function(parlist, m, n, q = 0)
  # Computes a_n for testing H_0 of m components
  # against H_1 of m+1 components
{
  if (q != 0) # an when the dimension of X is not zero.
    return (switch(as.character(q), "1" = 0.3, "2" = 2.0, "3" = 2.4, "4" = 2.4, 0.3))
  
  if (m == 1) {
    an <- 0.30
  }
  else if (m == 2) {
    omega <- omega.12(parlist)
    omega <- pmin(pmax(omega, 1e-16), 1-1e-16)  # an becomes NaN if omega[j]=0 or 1
    omega.term <- log(omega /(1-omega)) 
    
    # coefficients of -(intercept, misclterm, nterm, -atermcoeff^2)/atermcoeff
    b <- c(-4.937477, -0.845460, -56.494216, -0.21091) 
    x <- exp(b[1] + b[2] * omega.term + b[3] / n - log(2) / b[4])  # maxa=1
    an <- 0.25 * x / (1 + x)
    #   x <- exp(-1.642 - 0.434 * log(omega / (1 - omega)) - 101.80/n)  # maxa=2
    #   an <- 1.8 * x / (1 + x)
  }
  else if (m == 3) {
    omega <- omega.123(parlist)
    omega <- pmin(pmax(omega, 1e-16), 1-1e-16)  # an becomes NaN if omega[j]=0 or 1
    omega.12 <- omega[1]
    omega.23 <- omega[2]
    omega.term <- log(omega.12 * omega.23 / ((1-omega.12)*(1-omega.23)))

    b <- c(-2.4481555, -0.2034425, -56.9171687, -0.27104) 
    x <- exp(b[1] + b[2] * omega.term + b[3] / n - log(2) / b[4])  # maxa=1
    an <- 0.25 * x / (1 + x)
    # an <- 0.80 * x / (1 + x)
    #   x <- exp(-1.678 - 0.232 * log(t_omega) - 175.50/n)
    #   an <- 1.5 * x / (1 + x)
  } else if (m == 4) {
    omega <- omega.1234(parlist)
    omega <- pmin(pmax(omega, 1e-16), 1-1e-16)  # an becomes NaN if omega[j]=0 or 1
    omega.12 <- omega[1]
    omega.23 <- omega[2]
    omega.34 <- omega[3]
    omega.term <- log(omega.12 * omega.23 * omega.34 / 
                  ((1-omega.12)*(1-omega.23)*(1-omega.34)))
    b <- c(-5.3663749, -0.2462147, -199.6375112, -0.300460) 
    x <- exp(b[1] + b[2] * omega.term + b[3] / n - log(2) / b[4])  # maxa=1
    an <- 0.25 * x / (1 + x)
  }
  else 
    an <- 1.0
  
  return (an)
}  # end function anFormula
