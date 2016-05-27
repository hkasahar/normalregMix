\name{print.normalregMix}
\alias{print.normalregMix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prints normalregmix%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.normalregMix(x, ...)
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
    cat("\nCall:\n")
    print(x$call)
    callname <- deparse(x$call[1])
    if (grepl("MEMtest", callname)) {
        cat(sprintf("\nTesting the null hypothesis of \%d components\n", 
            x$m))
        cat("                            k = 1  k = 2  k = 3 \n")
        cat(c("modified EM-test statistic ", sprintf("\%.3f ", 
            x$emstat)), "\n")
        if (x$crit.method == "asy") {
            cat(c("asymptotic p-value         ", sprintf("\%.3f ", 
                x$pvals)), "\n")
        }
        else if (x$crit.method == "boot") {
            cat(c("bootstrap p-value          ", sprintf("\%.3f ", 
                x$pvals)), "\n")
        }
    }
    else if (grepl("PMLE", callname)) {
        cat("\nCoefficients:\n")
        print(x$coefficients, digits = 4)
    }
    else {
        stop("The object is not a valid normalMix object.")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
