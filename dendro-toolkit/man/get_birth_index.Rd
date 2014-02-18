\name{get_birth_index}
\alias{get_birth_index}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_birth_index(birth_year, tra, sparse)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{birth_year}{
%%     ~~Describe \code{birth_year} here~~
}
  \item{tra}{
%%     ~~Describe \code{tra} here~~
}
  \item{sparse}{
%%     ~~Describe \code{sparse} here~~
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
function (birth_year, tra, sparse) 
{
    if (sparse) {
        years <- sort(as.numeric(levels(tra$t)))
    }
    else {
        years <- sort(as.numeric(dimnames(tra)[[2]]))
    }
    year_index <- which(birth_year == years)
    age_index <- 1
    birth_index <- year_index - age_index
    return(birth_index)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
