\name{normalmixPMLE}
\alias{normalmixPMLE}
\title{
Penalized Maximum Likelihood Estimation of a Finite Mixture of Univariate Normals
}
\description{
Estimate parameters of a finite mixture of univariate normals by the method of penalized maximum likelhood.
}
\usage{
normalmixPMLE(y, m = 2, z = NULL, vcov.method = c("Hessian", "OPG", "none"), 
              ninits = 25, epsilon = 1e-08, maxit = 2000,
              epsilon.short = 0.01, maxit.short = 500)
}
\arguments{
  \item{y}{
  An n-vector of response values.
  }
  \item{m}{
  The number of components in the mixture.
  }
  \item{z}{
  An nxp matrix of covariates whose coefficients are common to all the components.
  }
  \item{vcov.method}{
  The method to compute the covariance matrix of the PMLE
  }
  \item{ninits}{
  The number of randomly drawn initial values.
  }
  \item{epsilon}{
  The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}. 
}
  \item{maxit}{
  The maximum number of iterations.
  }
  \item{epsilon.short}{
  The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
  }
  \item{maxit.short}{
  The maximum number of iterations in short EM.
  }
}
\details{
The likelihood function of normal mixtures diverges to infinity when the variance of a component is allowed to go to zero. The penalized MLE overcomes this problem by the use of a penalty term. See Chen et al. (2008) for details.

When \code{z=NULL}, \code{normalmixPMLE} estimates an \eqn{m}-component normal mixture model in which the \eqn{j}th component has mean \eqn{\mu_j} and variance \eqn{\sigma_j^2}. When \code{z} is not \code{NULL}, \code{normalmixPMLE} estimates an model in which the \eqn{j}th component has mean \eqn{\mu_j + \gamma'z} and variance \eqn{\sigma_j^2}. Here, \eqn{z} is a vector of covariates whose coefficient \eqn{\gamma} is common to all the components.

%%\deqn{
%%\sum_{j=1}^m \alpha_j N(\mu_j + \gamma z, \sigma_j^2), 
%%}
%%where \eqn{\alpha_j} is the mixing proportion of the \eqn{j}th component, and \eqn{f(y|z;\mu_j,\sigma_j^2,\gamma)} is the normal density with mean \eqn{\mu_j + \gamma z} and variance \eqn{\sigma_j^2}.
%%  ~~ If necessary, more details than the description above ~~


 }
\note{\code{normalmixPMLE} maximizes the penalized log-likelihood function using the EM algorithm with combining short and long runs of EM steps as in Biernacki et al. (2003). \code{normalmixPMLE} first runs the EM algorithm from \code{ninits}\eqn{* 4m(1 + p)} initial values with the convertence criterion \code{epsilon.short} and \code{maxit.short}. Then, \code{normalmixPMLE} uses \code{ninits} best initial values to run the EM algorithm with the convertence criterion \code{epsilon} and \code{maxit}.}
\value{
\code{normalmixPMLE} returns a list of class \code{normareglMix} with items:
\item{coefficients}{A vector of parameter estimates. Ordered as \eqn{\alpha_1,\ldots,\alpha_m,\mu_1,\ldots,\mu_m,\sigma_1,\ldots,\sigma_m,\gamma}.}
\item{parlist}{The parameter estimates as a list containing alpha, mu, and sigma (and gamma if z is included in the model).}
\item{vcov}{The estimated variance-covariance matrix.}
\item{loglik}{The maximized value of the log-likelihood.}
\item{penloglik}{The maximized value of the penalized log-likelihood.}
\item{aic}{Akaike Information Criterion of the fitted model.}
\item{bic}{Bayesian Information Criterion of the fitted model.}
\item{postprobs}{An nxm matrix of posterior probabilities for observations.}
\item{call}{The matched call.}
\item{m}{The number of components in the mixture.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  Biernacki, C., Celeux, G. and Govaert, G. (2003)
  Choosing Starting Values for the EM Algorithm for Getting the
          Highest Likelihood in Multivariate Gaussian Mixture Models,
  \emph{Computational Statistics and Data Analysis}, \bold{41}, 561--575.

  Chen, J., Tan, X. and Zhang, R. (2008)
  Inference for Normal Mixtures in Mean and Variance,
  \emph{Statistica Sinica}, \bold{18}, 443--465.
  
  McLachlan, G. J. and Peel, D. (2000) \emph{Finite Mixture Models}, John Wiley \& Sons, Inc.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--    or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
