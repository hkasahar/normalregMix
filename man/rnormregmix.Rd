\name{rnormregmix}
\alias{rnormregmix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Title %%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rnormregmix(n, x = NULL, alpha, mubeta, sigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{mubeta}{
%%     ~~Describe \code{mubeta} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
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
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (n, x = NULL, alpha, mubeta, sigma) 
{
    m <- length(alpha)
    mubeta <- matrix(mubeta, ncol = m)
    if (!is.null(x)) {
        x <- as.matrix(x)
        if (nrow(x) != n) {
            stop("y and x must have the same number of rows.")
        }
        x1 <- cbind(1, x)
        ii <- sample(m, n, replace = TRUE, prob = alpha)
        y <- rnorm(n, mean = rowSums(x1 * t(mubeta[, ii])), sd = sigma[ii])
    }
    else {
        ii <- sample(m, n, replace = TRUE, prob = alpha)
        y <- rnorm(n, mean = mubeta[, ii], sd = sigma[ii])
    }
    y
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
