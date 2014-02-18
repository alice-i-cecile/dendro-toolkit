\name{extract_cv_gnm}
\alias{extract_cv_gnm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extract_cv_gnm(growth_model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{growth_model}{
%%     ~~Describe \code{growth_model} here~~
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
function (growth_model) 
{
    coeff <- coef(growth_model)
    names(coeff) <- gsub("Mult.*)\\.", "", names(coeff))
    Q.cv <- coeff[grep("Q.cv", names(coeff))]
    F.cv <- coeff[grep("F.cv", names(coeff))]
    A.cv <- coeff[grep("A.cv", names(coeff))]
    names(Q.cv) <- substr(names(Q.cv), 5, nchar(names(Q.cv)))
    names(F.cv) <- substr(names(F.cv), 5, nchar(names(F.cv)))
    names(A.cv) <- substr(names(A.cv), 5, nchar(names(A.cv)))
    cv <- list(Q = Q.cv, F = F.cv, A = A.cv)
    cv <- lapply(cv, abs)
    return(cv)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
