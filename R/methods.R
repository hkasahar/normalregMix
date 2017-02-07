normalregMix.test.on <- FALSE
normalregMix.test.seed <- 8888577

#' @description Turns on/off the test mode.
#' @export
#' @title testMode
#' @name testMode
#' @description When the modified EM-algorithm is run, initial values are randomly created
#' based on the data given. If the test mode is turned on, these initial values
#' are going to be created with the random seed provided. This method would be suitable
#' for users who would like to replicate experiments. By default, the test mode is turned off.
#' @param on Option to turn on the test mode
#' @param seed The random seed to be used for initialization
#' @param hide.message Determines whether to print the current seed and status
testMode <- function(on = FALSE, seed = 8888577, hide.message = TRUE)
{
  unlockBinding("normalregMix.test.on", getNamespace("normalregMix"))
  unlockBinding("normalregMix.test.seed", getNamespace("normalregMix"))
  assign("normalregMix.test.on", on, getNamespace("normalregMix"))
  assign("normalregMix.test.seed", seed, getNamespace("normalregMix"))
  lockBinding("normalregMix.test.on", getNamespace("normalregMix"))
  lockBinding("normalregMix.test.seed", getNamespace("normalregMix"))

  if (!hide.message)
    print(paste("The test mode is currently",
                switch(as.character(normalregMix.test.on), "TRUE" = "ON", "FALSE" = "OFF"),
                "with seed",
                as.character(normalregMix.test.seed)))
}

#' @description Prints the multiple plots for diagnosis. 
#' @export
#' @title plotDiag
#' @name plotDiag
#' @param components n vector of the indices of components for each observation
#' @param y n vector of data that represents dependent variables
#' @param x n by q matrix of data that represent(s) covariates where q >= 1
#' @param m The number of components
plotDiag <- function(components, y = y, x = x, m = 2)
{
  dimx <- dim(as.matrix(x))[2]
  if ((dimx <= 1) || (is.null(dimx)))
    return (NULL)

  pivot.names <- as.character(seq(1, dimx))
  if (!is.null(colnames(x)))
    pivot.names <- colnames(x)

  ivs <- as.matrix(x)

  for (j in 1:m)
  {
    ivs.component <- as.matrix(ivs[components == j,])
    ys.component <- y[components == j]
    for (pivot in 1:dimx)
    {
      pivot.name <- pivot.names[pivot]
      ivs.pivot <- ivs.component[,pivot]
      ivs.others <- ivs.component[,-pivot]
      lm.y.other <- lm(ys.component ~ ivs.others)
      lm.pivot.other <- lm(ivs.pivot ~ ivs.others)
      plot.df <- data.frame(y.on.others = lm.y.other$residuals,
                            pivot.on.others = lm.pivot.other$residuals)
      plot <- ggplot(plot.df, aes(x=pivot.on.others, y=y.on.others))
      plot <- plot + geom_point(shape=1) + geom_smooth(method=lm) +
        xlab(paste("pivot.on.others (pivot on ", pivot.name, ", component ",
                   as.character(j), ")", sep = ""))
      print(plot)
    }
  }
}

#' @description Generates a vector that indicates which component each observation belongs to,
#' based on its posterior probability
#' @export
#' @title getComponentcomponents
#' @name getComponentcomponents
#' @param postprobs n by m matrix of posterior probabilities for
#' m-component model on n observations
#' @return n by 1 vector of components that indicate which component each observation belongs to
#' based on its posterior probability
getComponentcomponents <- function(postprobs)
{
  postprobs.mat <- as.matrix(postprobs)
  apply(postprobs.mat, 1, function(i) (which(i==max(i))))
}

#' @description Returns the summary of a normalregMix instance
#' @export
#' @title summary.normalregMix
#' @name summary.normalregMix
#' @param object normalregMix instance
#' @param reorder Determines whether components are reordered in summary
#' @param digits Digits used for reporting
#' @param ... Other arguments that do not affect the function
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

#' @description Prints the summary of a normalregMix instance
#' @export
#' @title print.summary.normalregMix
#' @name print.summary.normalregMix
#' @param x normalregMix instance
#' @param ... Other arguments that do not affect the function
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

#' @description Prints the description of a normalregMix instance
#' @export
#' @title print.normalregMix
#' @name print.normalregMix
#' @param x normalregMix instance
#' @param ... Other arguments that do not affect the function
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


#' @description Computes a_n based on empirical results found in Kasahara and Shimotsu (2015)
#' @export
#' @title anFormula
#' @name anFormula
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gamma_1, ..., gamma_m))
#' @param m The number of components in the mixture
#' @param n The number of observations
#' @param q The dimension of x (by default, 0)
#' @param LRT.penalized Determines whether penalized likelihood is used in calculation of LRT
#' statistic for likelihood in an alternative hypothesis.
#' @return a_n used to initialize values
#' @references Kasahara, H., and Shimotsu, K. (2015)
#' Testing the Number of Components in Normal Mixture Regression Models,
#' \emph{Journal of the American Statistical Association},
#' \bold{110}, 1632--1645.
anFormula <- function(parlist, m, n, q = 0, LRT.penalized = FALSE)
  # Computes a_n for testing H_0 of m components
  # against H_1 of m+1 components
{
  if (!LRT.penalized)
    return (anFormulaNotPenalized(parlist = parlist, m = m, n = n, q = q))
  if (q != 0) # an when the dimension of X is not zero.
    return (switch(as.character(q), "1" = 0.5, "2" = 2.0, "3" = 2.4, "4" = 2.4, 0.5))
  
  if (m == 1) {
    #   an <- 1.0
    an <- 0.50
    #   an <- 0.25
  }
  else if (m == 2) {
    omega <- omega.12(parlist)
    omega <- pmin(pmax(omega, 1e-16), 1-1e-16)  # an becomes NaN if omega[j]=0 or 1
    omega.term <- log(omega /(1-omega)) 
    
    # coefficients of -(intercept, misclterm, nterm, -atermcoeff^2)/atermcoeff
    b <- c(-3.3905358, -0.5152901, -39.1238054, -0.22429) 
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

#' @description Computes a_n based on empirical results found in Kasahara and Shimotsu (2015)
#' for the case where penalty term is not included in computation of LRT statistic.
#' @title anFormulaNotPenalized
#' @name anFormulaNotPenalized
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gamma_1, ..., gamma_m))
#' @param m The number of components in the mixture
#' @param n The number of observations
#' @param q The dimension of x (by default, 0)
#' @return a_n used to initialize values
#' @references Kasahara, H., and Shimotsu, K. (2015)
#' Testing the Number of Components in Normal Mixture Regression Models,
#' \emph{Journal of the American Statistical Association},
#' \bold{110}, 1632--1645.
anFormulaNotPenalized <- function(parlist, m, n, q = 0, ..)
{
  if (q != 0) # an when the dimension of X is not zero.
    return (switch(as.character(q), "1" = 0.01, "2" = 0.24, "3" = 0.31, "4" = 0.54, 0.5))
  
  if (m == 1) {
    an <- 0.50
  }
  else if (m == 2) {
    omega <- omega.12(parlist)
    omega <- pmin(pmax(omega, 1e-16), 1-1e-16)  # an becomes NaN if omega[j]=0 or 1
    omega.term <- log(omega /(1-omega)) 
    b <- c(-3.111849, -1.056635, -373.097585, -0.27019) # coefficients of -(intercept, misclterm, nterm, -atermcoeff^2)/atermcoeff
    x <- exp(b[1] + b[2] * omega.term + b[3] / n - log(2) / b[4])  # maxa=1
    an <- 1.5 * x / (1 + x)
    #   x <- exp(-1.642 - 0.434 * log(omega / (1 - omega)) - 101.80/n)  # maxa=2
    #   an <- 1.8 * x / (1 + x)
  }
  else if (m == 3) {
    omega <- omega.123(parlist)
    omega <- pmin(pmax(omega, 1e-16), 1-1e-16)  # an becomes NaN if omega[j]=0 or 1
    omega.12 <- omega[1]
    omega.23 <- omega[2]
    omega.term <- log(omega.12 * omega.23 / ((1-omega.12)*(1-omega.23)))
    
    b <- c(-1.3394916, -0.2054418, -191.9393302, -0.5673) # coefficients of -(intercept, misclterm, nterm, -atermcoeff^2)/atermcoeff
    x <- exp(b[1] + b[2] * omega.term + b[3] / n - log(2) / b[4])  # maxa=1
    an <- 1.5 * x / (1 + x)
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
    b <- c(-3.8800864, -0.2414154, -543.5967096, -0.4297) # coefficients of -(intercept, misclterm, nterm, -atermcoeff^2)/atermcoeff
    x <- exp(b[1] + b[2] * omega.term + b[3] / n - log(2) / b[4])  # maxa=1
    an <- 1.5 * x / (1 + x)
  }
  else 
    an <- 1.0
  
  return (an)
}  