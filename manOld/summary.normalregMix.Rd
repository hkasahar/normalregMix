\name{summary.normalregMix}
\alias{summary.normalregMix}
\title{Summarizing normal mixture and regression mixture model fits}
\description{
\code{\link[base]{summary}} method for class \code{normalregMix}.
}
\usage{
\method{summary}{normalregMix}(object, reorder = FALSE, digits = 4, \dots)
}
\arguments{
\item{object}{an object of class \code{normalregMix} such as a result of a call
to \code{\link{normalmixPMLE}}}
\item{reorder}{If FALSE, then the coefficients are ordered as \code{alpha}, \code{mu} or \code{mubeta}, \code{sigma}, \code{gamma}. If TRUE, then the coefficients are ordered component-by-component.
}
\item{digits}{Significant digits for printing values}
\item{\dots}{further arguments passed to \code{print} method.}
}
\details{
 \code{\link{summary.normalregMix}} prints parameter estimates, standard errors, p-values, etc. of a fitted mixture model.
}
\value{
The function \code{\link{summary.normalregMix}} returns a list of type \code{summary.normalregMix} with the following components:
\item{coefficients}{A vector of parameter estimates, reordered if \code{reorder}=TRUE.}
\item{parlist}{The parameter estimates as a list containing \code{alpha}, \code{mu} or \code{mubeta}, and \code{sigma} (and \code{gamma} if \code{z} is included in the model).}
\item{vcov}{The estimated variance-covariance matrix, reordered if \code{reorder}=TRUE.}
\item{loglik}{The maximized value of the log-likelihood.}
\item{penloglik}{The maximized value of the penalized log-likelihood.}
\item{aic}{Akaike Information Criterion of the fitted model.}
\item{bic}{Bayesian Information Criterion of the fitted model.}
\item{call}{The matched call.}
\item{m}{The number of components in the mixture.}
\item{reorder}{The value of \code{reorder}.}
\item{digits}{The value of \code{digits}.}
}
\seealso{
\code{\link{normalmixPMLE}},
\code{\link{regmixPMLE}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(faithful)
out <- normalmixPMLE(faithful$waiting, m=2)
summary(out)
summary(out, reorder=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
