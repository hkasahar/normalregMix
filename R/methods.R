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
    if (!is.null(object$parlist$gamma)) {
      p <- length(object$parlist$gamma)
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
  if (!is.null(x$parlist$gamma)) {
    p <- length(x$parlist$gamma)
    gamma <- tab[(nrow(tab)-p+1):nrow(tab), , drop = FALSE]
    cat("Gamma\n")
    printCoefmat(gamma, P.values = TRUE, has.Pvalue = TRUE, signif.legend = FALSE, digits=x$digits)
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

print.summary.normalregMix.old <- function(x, ...) {

m   <- x$m
tab <- x$coefficients
reorder <- x$reorder
k <- nrow(x$parlist$mubeta)
if (is.null(k)) { k <- 1 }

if (reorder) {
  cat(sprintf("Number of components: %d\n",m))
  for (j in 1:m){
    cat(sprintf("\nComponent %i\n",j))
    coef.j = tab[c(j, (m+(j-1)*k+1):(m+j*k), j+m*(k+1)), ]
    # rownames(coef.j) <- c("alpha","mu","sigma")
    printCoefmat(coef.j, P.values = TRUE, has.Pvalue = TRUE, signif.legend = FALSE, digits=x$digits)
  }
  if (!is.null(x$parlist$gamma)) {
    # browser()
    p <- length(x$parlist$gamma)
    gamma <- tab[(nrow(tab)-p+1):nrow(tab), , drop = FALSE]
    cat("\nGamma\n")
    printCoefmat(gamma, P.values = TRUE, has.Pvalue = TRUE, signif.legend = FALSE, digits=x$digits)
  }
  cat("---\n Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
} else {
  printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE, digits = x$digits)
  cat(sprintf("\nloglik at estimate: %.3f\n", x$loglik))
  cat(sprintf("penloglik at estimate: %.3f\n", x$penloglik))
  cat(sprintf("AIC: %.3f\n", x$aic))
  cat(sprintf("BIC: %.3f\n", x$bic))
}

}


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