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
H[,,1] <- Z
H[,,2] <- t(t(Z^2-1)/2/sigma^2)
H[,,3] <- t(t(Z^3-3*Z)/6/sigma^3)
H[,,4] <- t(t(Z^4-6*Z^2+3)/24/sigma^4)
return(H)
}  # end function hermite


tKR <- function (a,b) {
  # Computes the transpose of the Khatri-Rao product of t(a) and t(b)
  n <- nrow(a)
  k <- ncol(a)
  KR <- matrix(unlist(lapply(1:k, function(i) (a[,i]*b))),nrow=n)
  KR
}


coef.to.list <- function(coefficients, z = NULL) {
# 　Convert coefficients to list
len     <- length(coefficients)
p       <- 0
gamma   <- NULL

if (!is.null(z)) {
  z <- as.matrix(z)
  p <- ncol(z)
  gamma <- coefficients[(len-p+1):len]
}

m <- (len-p)/3
if (round(m) != m) {
  stop("The dimension of the coefficients is incompatible with z. Please check the data.")
}

param   <- matrix(coefficients[1:(len-p)], nrow=m, ncol=3)
alpha   <- param[, 1]
mu      <- param[, 2]
sigma   <- param[, 3]

a = list(alpha = alpha, mu = mu, sigma = sigma, gamma = gamma)

a

}
