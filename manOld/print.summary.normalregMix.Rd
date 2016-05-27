\name{print.summary.normalregMix}
\alias{print.summary.normalregMix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
print summary normalregmix%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.summary.normalregMix(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (x, ...) 
{
    m <- x$m
    tab <- x$coefficients
    reorder <- x$reorder
    k <- nrow(x$parlist$mubeta)
    if (is.null(k)) {
        k <- 1
    }
    if (reorder) {
        cat(sprintf("Number of components: \%d\n", m))
        for (j in 1:m) {
            cat(sprintf("\nComponent \%i\n", j))
            coef.j = tab[c(j, (m + (j - 1) * k + 1):(m + j * 
                k), j + m * (k + 1)), ]
            printCoefmat(coef.j, P.values = TRUE, has.Pvalue = TRUE, 
                signif.legend = FALSE)
        }
        if (!is.null(x$parlist$gamma)) {
            p <- length(x$parlist$gamma)
            gamma <- tab[(nrow(tab) - p + 1):nrow(tab), , drop = FALSE]
            cat("\nGamma\n")
            printCoefmat(gamma, P.values = TRUE, has.Pvalue = TRUE, 
                signif.legend = FALSE)
        }
        cat("---\n Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
    }
    else {
        printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE, 
            digits = 4)
        cat(sprintf("\nloglik at estimate: \%.3f\n", x$loglik))
        cat(sprintf("penloglik at estimate: \%.3f\n", x$penloglik))
        cat(sprintf("AIC: \%.3f\n", x$aic))
        cat(sprintf("BIC: \%.3f\n", x$bic))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
