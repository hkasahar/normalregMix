#' @description Computes the normalized Hermite polynomial for computing
#' the critical values and p-values of the modified EM test.
#' @export
#' @title hermite
#' @name hermite
#' @param Z n by m matrix of normalized data.
#' @param sigma m by 1 vector of parameters.
#' @return n by m by 4 array of normalized Hermite polynomials.
#' @keywords{intenal}
hermite <- function(Z, sigma)
# Computes the normalized Hermite polynomial for computing
# the critical values and p-values of the modified EM test
# Input
#   Z (n by m) : normalized data
#  sigma (m by 1): parameters
# Output
#   H (n by m by 4): normalized Hermite polynomials
{
n <- nrow(Z)
m    <- length(sigma)

H <- array(0,dim=c(n,m,4))
H[,,1] <- Z/sigma
H[,,2] <- t(t(Z^2-1)/2/sigma^2)
H[,,3] <- t(t(Z^3-3*Z)/6/sigma^3)
H[,,4] <- t(t(Z^4-6*Z^2+3)/24/sigma^4)
return(H)
}  # end function hermite

#' @description Computes the transpose of the Khatri-Rao product of \eqn{a'} and \eqn{b'}.
#' @export
#' @title tKR
#' @name tKR
#' @param a n by k matrix.
#' @param b n by k matrix.
#' @return n by k*k matrix of the transpose of the Khatri-Rao product 
#' of \eqn{a'} and \eqn{b'}.
#' @keywords{intenal}
tKR <- function (a,b) {
  # Computes the transpose of the Khatri-Rao product of t(a) and t(b)
  n <- nrow(a)
  k <- ncol(a)
  KR <- matrix(unlist(lapply(1:k, function(i) (a[,i]*b))),nrow=n)
  KR
}

