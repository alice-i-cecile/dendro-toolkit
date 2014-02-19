\name{make_gam_formula}
\alias{make_gam_formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make_gam_formula(incQ, incF, incA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{incQ}{
%%     ~~Describe \code{incQ} here~~
}
  \item{incF}{
%%     ~~Describe \code{incF} here~~
}
  \item{incA}{
%%     ~~Describe \code{incA} here~~
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
function (incQ, incF, incA) 
{
    dep.str <- "G"
    ind.str <- "0"
    if (incQ) {
        ind.str <- paste(ind.str, "Q.cv", sep = "+")
    }
    if (incF) {
        ind.str <- paste(ind.str, "F.cv", sep = "+")
    }
    if (incA) {
        ind.str <- paste(ind.str, "lo(A.num, span=span, degree=degree)", 
            sep = "+")
    }
    formula.str <- paste(dep.str, ind.str, sep = "~")
    return(formula.str)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
